from GenomicFeatures import cer_db, alb_db
import SummaryTable
import glob
import os
import csv
import cPickle
import math
import pandas as pd
import itertools
import goatools, goatools.associations
import re
import Shared

from sklearn.cross_validation import StratifiedKFold
from sklearn.metrics import roc_curve, auc
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier

import numpy as np 
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3

# Cache the essentials for futher use:
CER_ESSENTIALS, CER_NON_ESSENTIALS = SummaryTable.get_cerevisiae_essentials()

def get_hits_from_bed(bed_file):
    """Read a BED-formatted hit file into a data structure acceptable by
    `SummaryTable.analyze_hits`.
    
    Returns
    -------
    list of dicts
        A list of hit objects, as expected by `SummaryTable.analyze_hits`. 
    """
    
    result = []
    with open(bed_file, 'r') as in_file:
        in_file.readline() # Drop the header header
        for line in in_file:
            chrom_name, start, _stop, _strand, score = line.split() 
            if chrom_name == 'chrM':
                continue
            hit_pos = int(start)
            
            fs = cer_db[chrom_name][hit_pos]
            if len(fs) == 0:
                name = ig_type = "nan"
            else:
                f = fs[0]
                name = f.standard_name
                ig_type = "ORF"
            
            result.append({"chrom": chrom_name, "hit_pos": hit_pos, "gene_name": name,
                           "ig_type": ig_type, "hit_count": (int(score)-100)/20})
            
    return result

def scale_cf(width=1, height=1,figure=None):
    if figure is None:
        figure = plt.gcf()
    if width != 1:
        figure.set_figwidth(plt.gcf().get_figwidth()*width)
    if height != 1:
        figure.set_figheight(plt.gcf().get_figheight()*height)


def test_classifier(classifier, features, annotations):
    """Test a classifier using given features and annotations.
    
    Parameters
    ----------
    classifier : sklearn Classifier
        The classifier to test.
    features : list of sequences
        A list of features per gene.
    annotations : list of 1 or 0
        A list of annotations for the genes.
    
    Returns
    -------
    tuple
        The ROC curve parameters, as returned by `sklearn.metrics.roc_curve`.
    """
    
    scores = np.empty(len(annotations))
    skf = StratifiedKFold(y=annotations, n_folds=5, shuffle=True, random_state=0)
    for roc_train, roc_test in skf:
        roc_train_set = [features[i] for i in roc_train]
        roc_test_set = [features[i] for i in roc_test]
        roc_train_essentiality = [annotations[i] for i in roc_train]
        
        classifier.fit(roc_train_set, roc_train_essentiality)
        scores[roc_test] = classifier.predict_proba(roc_test_set)[:,1]
    
    return roc_curve(annotations, scores)

def train_classifiers(training_records, essentials, feature_groups, classifier_factory):
    """Train given classifiers using given records and feature groups.
    
    Parameters
    ----------
    training_records : list of dicts
        A list of record objects, from which the features will be taken.
    essentials : set
        A set of standard gene names which are considered essential within the
        given record list. All records not in this set are considered
        non-essential. 
    feature_groups : sequence of (str, sequence)
        A sequence of all feature groups to train with. Every feature group
        should be represented by a tuple of
        (feature group name, [feature name sequence]).
    classifier_factory : dict of str to callable
        A dict mapping classifier names to their creator callables.
        
    Returns
    -------
    dict of tuple to tuple
        A mapping between (classifier class name, feature group name) to
        (classifier object, roc_curve_stats), where the roc_curve_stats are in
        the format returned by `sklearn.metrics.roc_curve`.
    """
    
    result = {}
    
    training_essentials = [1 if r["feature"].standard_name in essentials else 0
                           for r in training_records]
    
    for group_name, features in feature_groups:
        training_features = [[r[f] for f in features] for r in training_records]
        for cls_name, cls_creator in classifier_factory.items():
            # Fit the data in parts to get the ROC curve:
            fpr, tpr, thresholds = test_classifier(cls_creator(), training_features, training_essentials)
            
            # Fit the entirety of the data for the classification later on
            classifier = cls_creator()
            classifier.fit(training_features, training_essentials)
            result[(cls_name, group_name)] = (classifier, (fpr, tpr, thresholds))
            
    return result

def predict_records(records_to_predict, feature_groups, classifiers):
    """Predict essentiality in given records using the give feature groups and
    pre-trained classifiers.
    
    This function does not return anything, rather it injects the predicted
    values into the records.
    
    Parameters
    ----------
    records_to_predict : list of dict
        A list of record objects to perform the predictions on.
    feature_groups : sequence of (str, sequence)
        A sequence of all feature groups to train with. Every feature group
        should be represented by a tuple of
        (feature group name, [feature name sequence]).
    classifiers : dict of tuple to tuple
        The trained classifiers, in the format returned by `train_classifiers`.
    """
    
    for group_name, features in feature_groups:
        prediction_features = [[r[f] for f in features] for r in records_to_predict]
        for (cls_name, cls_group_name, _rdf), (classifier, _stats) in classifiers.items():
            if cls_group_name != group_name:
                continue
            predictions = classifier.predict_proba(prediction_features)
            pred_field = "%s-%s" % (cls_name, group_name)
            for record, prediction in zip(records_to_predict, predictions[:,1]):
                record[pred_field] = prediction

def get_training_groups():
    """Get all interesting permutations of features to train on.
    
    Note
    ----
    
    Use with caution, as this can be quite large.
    
    Returns
    -------
    generator of tuples
        Each tuple contains a permutation of features to use for training.
    """
    
    subgroups = (("neighborhood_index", "insertion_index"),
                 ("length", "hits", ("length", "hits",), ""),
                 ("max_free_region", "freedom_index", "adj_max_free_region", "adj_freedom_index", ""),
                 ("upstream_hits_50", "upstream_hits_100", ("upstream_hits_50", "upstream_hits_100"), ""),
                 ("std_dev", "index_of_dispersion", ""))
    
    def _flatten(features):
        # Recursive, generator-based flattening.
        for f in features:
            if isinstance(f, str):
                yield f
            else:
                for sub_f in f:
                    yield sub_f
    
    for group in itertools.product(*subgroups):
        yield tuple(f for f in _flatten(group) if f)

def test_training_data(records, essentials, output_folder):
    """Test the training records and store the results in an output folder.
    
    The output includes the FPR/TRP/threshold table, a table with the feature
    values and the predictions, the ROC curve plots, and general stats about
    the different classifiers.
    
    Parameters
    ----------
    records : list of dict
        A sequence of record objects.
    essentials : set of str
        A set of essential gene names (in their standard form). Genes not
        present in this set are considered non-essential.
    output_folder : str
        The path to the folder in which to store the output
    """
    
    Shared.make_dir(output_folder)
    
    classifier_factory = {
        "LR": lambda: LogisticRegression(), #(class_weight='auto')
        "RF": lambda: RandomForestClassifier(n_estimators=100, random_state=0)
    }
    
    stats_csv = csv.writer(open(os.path.join(output_folder, "stats.csv"), 'w'), delimiter=',')
    stats_csv.writerow(["Features"] + sum([["%s - %s" % (cname, param) for param in ("AUC", "Youden - TP", "Youden - FP", "Youden - threshold")] for cname in classifier_factory], []))
    
    # You can choose between manual control of the feature groups and the
    # auto-generated permutations:
    for columns_for_classification in (
                                        ("neighborhood_index", "length", "hits", "freedom_index", "upstream_hits_100"),
                                      ):
#     for columns_for_classification in get_training_groups():
        all_data_features = [[record[col_name]  for col_name in columns_for_classification] for record in records]
        all_essentiality_list = [1 if record["feature"].standard_name in essentials else 0 for record in records]
        
        folds = 5
        skf = StratifiedKFold(y=all_essentiality_list, n_folds=folds, shuffle=True, random_state=0)
        
        scale_cf(width=2)
        
        csv_row = []
        
        for ix, (classifier_name, classifier_creator) in enumerate(classifier_factory.iteritems()):
            scores = np.empty(len(records))
            
            for train, test in skf:
                classifier = classifier_creator()
                
                train_set = [all_data_features[i] for i in train]
                test_set = [all_data_features[i] for i in test]
                
                train_essentiality = [all_essentiality_list[i] for i in train]
                
                classifier.fit(train_set, train_essentiality)
                scores[test] = classifier.predict_proba(test_set)[:,1]
            
            plt.subplot(1, 2, ix+1) # TODO: will actually work with at most 2 classifiers!
            fpr, tpr, thresholds = roc_curve(all_essentiality_list, scores)
            classifier_auc = auc(fpr, tpr)
            plt.plot(fpr, tpr, label="%s (AUC:%.2f)" % (classifier_name, auc(fpr, tpr)))
            plt.legend(loc="lower right")
            
            y_threshold, y_tpr, y_fpr = get_youden_statistic(thresholds, tpr, fpr)
            csv_row.extend([classifier_auc, y_tpr, y_fpr, y_threshold])
            
            concat_features_str = "_".join(columns_for_classification)
            
            with open(os.path.join(output_folder, "%s.%s.csv" % (classifier_name, concat_features_str)), 'w') as roc_file:
                roc_writer = csv.writer(roc_file, delimiter=',')
                roc_writer.writerow(["Threshold", "TPR", "FPR"])
                for threshold, t, f in  zip(thresholds, tpr, fpr):
                    roc_writer.writerow([threshold, t, f])
                    
            with open(os.path.join(output_folder, "train_table.%s.%s.csv" % (classifier_name, concat_features_str)), 'w') as train_table:
                train_writer = csv.writer(train_table)
                train_writer.writerow(["Name"] + list(columns_for_classification) + ["Prediction", "Essential in literature"])
                for record, features, prediction, is_essential in zip(records, all_data_features, scores, all_essentiality_list):
                    train_writer.writerow([record["feature"].name] + features + [prediction, "Yes" if is_essential else "No"])
            
        stats_csv.writerow([", ".join(columns_for_classification)] + csv_row)
            
        plt.savefig(os.path.join(output_folder, 'ROCs.%s.pdf' % ('_'.join(columns_for_classification))),
                    transparent=True,
                    dpi=300)
        plt.close()

def get_youden_statistic(thresholds, tpr, fpr):
    # TODO: refactor in to the standard fpr/trp/threshold order.
    return max(zip(thresholds, tpr, fpr), key=lambda (th, t, f): t-f)

def classify_albicans(pre_prediction_hits, post_prediction_hits, train_hits,
                      train_essentials, output_folder, target_fpr=0.075,
                      run_go=False):
    """Classify albicans data.
    
    Parameters
    ----------
    pre_prediction_hits : list of dict
        A list of hits from the pre-evolution experiments.
    post_prediction_hits : list of dict
        A list of hits from the post-evolution experiments.
    train_hits : list of dict
        A list of the cerevisiae hits to train on.
    train_essentials : set of str
        A set of standard gene names that denote the essentials among the
        training dataset.
    output_folder : str
        The output folder which all results will be saved to.
    target_fpr : float
        The FPR which to choose when calling essentiality.
    run_go : bool
        A flag for running the GO analysis (can take some time).
    """
    
    # For convenience, if the output folder doesn't exist, create it:
    Shared.make_dir(output_folder)
    
    # TODO: ignore features with no hits in neighborhood!
    coverages = SummaryTable.get_alb_coverage()
    min_coverage = 0.95 # Min coverage required to be factored into the analysis
    # Use all names because we sometimes use the standard name and sometimes the CGDID:
    ignored_genes = set.union(*[f.all_names for f in alb_db.get_all_features() if coverages.get(f.standard_name, 1) < min_coverage])
    
    # Analyze the hits from the prediction datasets according to given RDFs:
    post_analyzed_datasets = {}
    pre_analyzed_datasets = {}
#     rdfs = [1, 5, 10, 15, 20, 25, 30] # Read depth filters
    rdfs = [1]
    for read_depth_filter in rdfs:
        filtered_dataset = SummaryTable.analyze_hits([o for o in post_prediction_hits if o["hit_count"] >= read_depth_filter], alb_db, 10000)
        post_analyzed_datasets[read_depth_filter] = [r for r in filtered_dataset.values() if r["feature"].is_orf and r["feature"].standard_name not in ignored_genes]
        
        filtered_dataset = SummaryTable.analyze_hits([o for o in pre_prediction_hits if o["hit_count"] >= read_depth_filter], alb_db, 10000)
        pre_analyzed_datasets[read_depth_filter] = [r for r in filtered_dataset.values() if r["feature"].is_orf and r["feature"].standard_name not in ignored_genes]
    
    # Define the classificaiton parameters:
    feature_groups = (
                     ("G3", ("neighborhood_index", "length", "hits", "freedom_index", "upstream_hits_100")),
                     )
    
    classifier_factory = {
        "LR": lambda: LogisticRegression(),
        "RF": lambda: RandomForestClassifier(n_estimators=100, random_state=0)
    }
    
    cls_names = classifier_factory.keys()
    group_names = [p[0] for p in feature_groups]
    cls_keys = [(cls_type, group) for cls_type in cls_names for group in group_names]
    
    # Train the classifiers:
    classifiers = {}
    for rdf in rdfs:
        train_records = SummaryTable.analyze_hits([h for h in train_hits if h["hit_count"] >= rdf], cer_db, 10000).itervalues()
        train_records = list( filter_cer_training_data(train_records) )
        rdf_classifiers = train_classifiers(train_records, train_essentials, feature_groups, classifier_factory)
        for (cls_name, grp_name), cls_data in rdf_classifiers.items():
            classifiers[(cls_name, grp_name, rdf)] = cls_data 
    
    # Call essentiality scores of the albicans datasets:
    for predict_datasets in (pre_analyzed_datasets, post_analyzed_datasets):
        for rdf, dataset in predict_datasets.items():
            rdf_classifiers = dict(item for item in classifiers.items() if item[0][2] == rdf)
            predict_records(dataset, feature_groups, rdf_classifiers)
            
    # Find the thresholds for essentiality calling:
    threshold_ixs = {}
    thresholds = {}
    for (cls_name, grp_name, rdf), (_cls, (fpr, tpr, ts)) in classifiers.iteritems():
        closest_ix = np.absolute(fpr - target_fpr).argmin()
        threshold_ixs[(cls_name, grp_name, rdf)] = closest_ix
        thresholds[(cls_name, grp_name, rdf)] = ts[closest_ix]

    # Inject some of the pre-evolution data into the post-evolution records for
    # easier display:
    for read_depth_filter, post_dataset in post_analyzed_datasets.iteritems():
        pre_dataset = pre_analyzed_datasets[read_depth_filter]
        for post_record in post_dataset:
            pre_record = (r for r in pre_dataset if r["feature"] == post_record["feature"]).next()
            
            for cls in cls_names:
                for grp, _cols in feature_groups:
                    key = "%s-%s" % (cls, grp)
                    pre_score = pre_record[key]
                    post_score = post_record[key]
                    
                    post_record["pre-" + key] = pre_score
                    
                    threshold = thresholds[(cls, grp, read_depth_filter)]
                    if pre_score >= threshold:
                        if post_score >= threshold:
                            verdict = "Essential"
                        else:
                            verdict = "Weirdo"
                    else:
                        if post_score >= threshold:
                            verdict = "Sick"
                        else:
                            verdict = "Non-essential"
                            
                    post_record[key + "-verdict"] = verdict
            
            # Log2 of NI is tricky because the NI can be zero. It may be
            # worthwhile to adjust it so it's never absolute zero, but right
            # now it's unclear how to do this. Meanwhile, we use arbitrary
            # values for the cases which we can't compute.
            pre_ni = pre_record["neighborhood_index"]
            post_ni = post_record["neighborhood_index"]
            post_record["pre_neighborhood_index"] = pre_ni
            if pre_ni == 0 and post_ni == 0:
                log2ni = 0
            elif pre_ni == 0:
                log2ni = 999
            elif post_ni == 0:
                log2ni = -999
            else:
                log2ni = math.log(post_ni / pre_ni, 2)
            post_record["log2_ni"] = log2ni
            
            post_record["s_coefficient"] = post_record["s_value"] - pre_record["s_value"]
            
    cols_config = [
        {
            "field_name": "feature",
            "csv_name": "Standard name",
            "format": lambda f: f.standard_name
        },
        
        {
            "field_name": "feature",
            "csv_name": "Common name",
            "format": lambda f: f.common_name
        },
        
        {
            "field_name": "feature",
            "csv_name": "Sc ortholog",
            "format": lambda f: ','.join(f.cerevisiae_orthologs)
        },
        
        {
            "field_name": "cer_fitness",
            "csv_name": "Sc fitness",
        },
        
        {
            "field_name": "essential_in_cerevisiae",
            "csv_name": "Essential in Sc",
        },
        
        {
            "field_name": "cer_synthetic_lethal",
            "csv_name": "Synthetic lethal in Sc",
        },
        
        {
            "field_name": "pombe_ortholog",
            "csv_name": "Sp ortholog",
        },
         
        {
            "field_name": "essential_in_pombe",
            "csv_name": "Essential in Sp",
        },
        
        {
            "field_name": "essential_in_albicans",
            "csv_name": "Essential in Calb",
        },
        
        {
            "field_name": "essential_in_albicans_grace_roemer",
            "csv_name": "Essential in Calb - GRACE - Roemer",
        },
        
        {
            "field_name": "essential_in_albicans_grace_omeara",
            "csv_name": "Essential in Calb - GRACE - O'Meara",
        },
        
        {
            "field_name": "hits",
            "csv_name": "Hits",
            "format": "%d"
        },
                   
        {
            "field_name": "feature",
            "csv_name": "Length",
            "format": lambda f: len(f)
        },
        
        {
            "field_name": "max_free_region",
            "csv_name": "Max free region",
            "format": "%d"
        },
                   
        {
            "field_name": "freedom_index",
            "csv_name": "Freedom index",
            "format": "%.2f"
        },
        
        {
            "field_name": "neighborhood_index",
            "csv_name": "Neighborhood index",
            "format": "%.3f"
        },
        
        {
            "field_name": "log2_ni",
            "csv_name": "Log2 NI (post/pre)",
            "format": "%.1f"
        },
        
        {
            "field_name": "log2_ni_adj",
            "csv_name": "Log2 NI adjusted",
            "format": "%.2f"
        },
                   
        {
            "field_name": "s_coefficient",
            "csv_name": "S coefficient (log2)",
            "format": "%.2f"
        },
        
        {
            "field_name": "upstream_hits_100",
            "csv_name": "Upstream hits 100",
            "format": "%d"
        },
        
        ] + \
        [{"field_name": "%s%s-%s" % (ds, cls, grp),
          "csv_name": "%s - %s - %s" % ({"pre-": "Pre", "": "Post"}[ds], cls, grp),
          "format": "%.3f"}
         for cls in cls_names for grp in group_names for ds in ("pre-", "")] + \
        [{"field_name": "%s-%s-verdict" % (cls, grp),
          "csv_name": "%s - %s - Verdict" % (cls, grp)}
         for cls in cls_names for grp in group_names
        ] + \
        [
        
        {
            "field_name": "unique_coverage",
            "csv_name": "Unique coverage",
            "format": "%.2f"
        },
        
        {
            "field_name": "feature",
            "csv_name": "Type",
            "format": lambda f: f.type
        },
        
        {
            "field_name": "feature",
            "csv_name": "Description",
            "format": lambda f: f.description
        },
    ]
    
    # TODO: the run_go flag needs to be refactored.  
    if run_go:
        # TODO: do we really have to use CGDID? Why can't use the standard name?
        # Initialize the GO-related data:
        obodag = goatools.obo_parser.GODag(Shared.get_dependency("albicans/go-basic.obo"))
        gene2go = goatools.associations.read_associations(Shared.get_dependency("albicans/gene2go.txt"))
        bg_orfs = set(f.primary_cgdid for f in alb_db.get_all_features() if f.is_orf and f.primary_cgdid not in ignored_genes)
        # GO analysis requires a set of background genes. It may make sense to
        # test against different backgrounds. This analysis object is constructed
        # with a background of all tested albicans ORFs.
        goeaobj_vs_all = goatools.go_enrichment.GOEnrichmentStudy(
            bg_orfs,
            gene2go, # geneid/GO associations
            obodag, # Ontologies
            propagate_counts = False,
            alpha = 0.05, # default significance cut-off
            methods = ['fdr_bh'] # defult multipletest correction method
        )
    else:
        goeaobj_vs_all = None
    
    def run_go_on_gene_list(gene_list, goeaobj, filename):
        cgdid_list = [alb_db.get_feature_by_name(f).primary_cgdid for f in gene_list]
        with open(os.path.join(output_folder, filename), 'w') as out_file:
            out_file.write("\n".join(cgdid_list))
        
        if not run_go:
            return
        
        base_filename = os.path.join(output_folder, os.path.splitext(filename)[0])
        goea_results = [r for r in goeaobj.run_study(cgdid_list) if r.p_fdr_bh < 0.05]
        goeaobj.wr_tsv(base_filename + ".go.tsv", goea_results)
    
    # Output all our analysis results in human-readable ways:
    for read_depth_filter, dataset in post_analyzed_datasets.iteritems():
        output_suffix = "filter_%d" % read_depth_filter
        output_filename = os.path.join(output_folder, "albicans_essentiality_prediction.%s.csv" % output_suffix)
        
        SummaryTable.enrich_alb_records(dataset)
        
        SummaryTable.write_data_to_csv(dataset, cols_config, output_filename)
        
        alb_orfs = set(r["feature"].standard_name for r in dataset)
        alb_literature_essentials = set(r["feature"].standard_name for r in dataset if r["essential_in_albicans"] == "Yes")
        alb_literature_non_essentials = set(r["feature"].standard_name for r in dataset if r["essential_in_albicans"] == "No")
        alb_omeara_essentials = set(r["feature"].standard_name for r in dataset if r["essential_in_albicans_grace_omeara"] == "Yes")
        alb_omeara_non_essentials = set(r["feature"].standard_name for r in dataset if r["essential_in_albicans_grace_omeara"] == "No")
        alb_omeara_orfs = alb_omeara_essentials | alb_omeara_non_essentials

        cer_orthologs = set(r["feature"].standard_name for r in dataset if r["feature"].cerevisiae_orthologs)
        pom_orthologs = set(r["feature"].standard_name for r in dataset if r["pombe_ortholog"])
        cer_essentials = set(r["feature"].standard_name for r in dataset if r["essential_in_cerevisiae"] == "Yes" or r["cer_synthetic_lethal"] == "Yes")
        pom_essentials = set(r["feature"].standard_name for r in dataset if r["essential_in_pombe"] == "Yes")

        ortholog_intersection = cer_orthologs & pom_orthologs
        
        if run_go:
            # A GO analysis object that considers the background to be all genes
            # with orthologs in cerevisiae and pombe:
            goeaobj_vs_orthologs = goatools.go_enrichment.GOEnrichmentStudy(
                set(f.primary_cgdid for f in [alb_db.get_feature_by_name(fn) for fn in ortholog_intersection]),
                gene2go, # geneid/GO associations
                obodag, # Ontologies
                propagate_counts = False,
                alpha = 0.05, # default significance cut-off
                methods = ['fdr_bh'] # defult multipletest correction method
            )
        else:
            goeaobj_vs_orthologs = None
        
        draw_venn("Orthologs",
                  os.path.join(output_folder, "orthologs_%d.png" % read_depth_filter),
                  [alb_orfs, cer_orthologs, pom_orthologs],
                  ("Calb", "Sc", "Sp"))
        
        draw_venn("Literature essentials",
                  os.path.join(output_folder, "lit_essentials_%d.png" % read_depth_filter),
                  [alb_literature_essentials, cer_essentials, pom_essentials],
                  ("Calb", "Sc", "Sp"))
        
        # Create the Venn diagrams:
        for classifier_type in cls_names:
            for group_name, _cols in feature_groups:
                threshold = thresholds[(classifier_type, group_name, read_depth_filter)]
                alb_essential_records = [r for r in dataset if r["%s-%s" % (classifier_type, group_name)] >= threshold and 
                                         r["pre-%s-%s" % (classifier_type, group_name)] >= threshold]
                alb_classifier_essentials = set(r["feature"].standard_name for r in alb_essential_records)
                alb_classifier_non_essentials = alb_orfs - alb_classifier_essentials
                
                # Write a comparison table vs. GRACE:
                with open(os.path.join(output_folder, "predicted_vs_grace.csv"), 'w') as out_file:
                    writer = csv.writer(out_file)
                    writer.writerow(["", "GRACE ess.", "GRACE non ess.", "Total"])
                    writer.writerow(["Predicted ess.",
                                     len(alb_classifier_essentials & alb_omeara_essentials),
                                     len(alb_classifier_essentials & alb_omeara_non_essentials),
                                     len(alb_classifier_essentials)])
                    writer.writerow(["Predicted non ess.",
                                     len(alb_classifier_non_essentials & alb_omeara_essentials),
                                     len(alb_classifier_non_essentials & alb_omeara_non_essentials),
                                     len(alb_classifier_non_essentials)])
                    writer.writerow(["Total",
                                     len(alb_omeara_essentials),
                                     len(alb_omeara_non_essentials)])
                
                draw_venn("Essentials (pred. vs. O'Meara) - %s, %s, %.2f, %d" % (classifier_type, group_name, threshold, read_depth_filter),
                          os.path.join(output_folder, "pred_vs_omeara_essentials_%s_%s_%.2f_%d.png" % (classifier_type, group_name, threshold, read_depth_filter)),
                          [alb_classifier_essentials & alb_omeara_orfs, alb_omeara_essentials],
                          ("Predicted", "O'Meara essentials"))
                
                draw_venn("Non-essentials (pred. vs. O'Meara) - %s, %s, %.2f, %d" % (classifier_type, group_name, threshold, read_depth_filter),
                          os.path.join(output_folder, "pred_vs_omeara_non_essentials_%s_%s_%.2f_%d.png" % (classifier_type, group_name, threshold, read_depth_filter)),
                          [alb_classifier_non_essentials & alb_omeara_orfs, alb_omeara_non_essentials],
                          ("Predicted non-essentials", "O'Meara non-essentials"))
                
                draw_venn("Essentials (pred. vs. lit.) - %s, %s, %.2f, %d" % (classifier_type, group_name, threshold, read_depth_filter),
                          os.path.join(output_folder, "pred_vs_lit_essentials_%s_%s_%.2f_%d.png" % (classifier_type, group_name, threshold, read_depth_filter)),
                          [alb_classifier_essentials, alb_literature_essentials],
                          ("Predicted", "Literature"))
                
                draw_venn("Essentials (pred. vs. lit. non-ess.) - %s, %s, %.2f, %d" % (classifier_type, group_name, threshold, read_depth_filter),
                          os.path.join(output_folder, "pred_vs_lit_non_essentials_%s_%s_%.2f_%d.png" % (classifier_type, group_name, threshold, read_depth_filter)),
                          [alb_classifier_non_essentials, alb_literature_non_essentials],
                          ("Predicted non-essentials", "Literature non-essentials"))
                
                alb_essentials_with_orthologs = alb_classifier_essentials & ortholog_intersection
                
                draw_venn("Non-essentials (orthologs) - %s, %s, %.2f, %d" % (classifier_type, group_name, threshold, read_depth_filter),
                          os.path.join(output_folder, "non_essentials_orthologs_%s_%s_%.2f_%d.png" % (classifier_type, group_name, threshold, read_depth_filter)),
                          [alb_classifier_non_essentials & ortholog_intersection,
                           cer_orthologs - cer_essentials,
                           pom_orthologs - pom_essentials],
                          ("Calb", "Sc", "Sp"))
                
                draw_venn("Essentials (orthologs) - %s, %s, %.2f, %d" % (classifier_type, group_name, threshold, read_depth_filter),
                          os.path.join(output_folder, "essentials_orthologs_%s_%s_%.2f_%d.png" % (classifier_type, group_name, threshold, read_depth_filter)),
                          [alb_essentials_with_orthologs, cer_essentials, pom_essentials],
                          ("Calb", "Sc", "Sp"))
                
                gene_list_filename_template = "%s_%s_%s_%%s.txt" % (classifier_type, group_name, read_depth_filter)
                run_go_on_gene_list(alb_essentials_with_orthologs - (cer_essentials | pom_essentials), goeaobj_vs_all, gene_list_filename_template % ("alb"))
                run_go_on_gene_list(cer_essentials - (alb_essentials_with_orthologs | pom_essentials), goeaobj_vs_all, gene_list_filename_template % ("cer"))
                run_go_on_gene_list(pom_essentials - (cer_essentials | alb_essentials_with_orthologs), goeaobj_vs_all, gene_list_filename_template % ("pom"))
                run_go_on_gene_list((alb_essentials_with_orthologs & cer_essentials) - pom_essentials, goeaobj_vs_all, gene_list_filename_template % ("alb_cer"))
                run_go_on_gene_list((alb_essentials_with_orthologs & pom_essentials) - cer_essentials, goeaobj_vs_all, gene_list_filename_template % ("alb_pom"))
                run_go_on_gene_list((cer_essentials & pom_essentials) - alb_essentials_with_orthologs, goeaobj_vs_all, gene_list_filename_template % ("cer_pom"))
                run_go_on_gene_list(alb_essentials_with_orthologs & cer_essentials & pom_essentials, goeaobj_vs_all, gene_list_filename_template % ("alb_cer_pom"))
                
                run_go_on_gene_list(alb_essentials_with_orthologs - (cer_essentials | pom_essentials), goeaobj_vs_orthologs, gene_list_filename_template % ("alb_vs_orth"))
                run_go_on_gene_list(cer_essentials - (alb_essentials_with_orthologs | pom_essentials), goeaobj_vs_orthologs, gene_list_filename_template % ("cer_vs_orth"))
                run_go_on_gene_list(pom_essentials - (cer_essentials | alb_essentials_with_orthologs), goeaobj_vs_orthologs, gene_list_filename_template % ("pom_vs_orth"))
                run_go_on_gene_list((alb_essentials_with_orthologs & cer_essentials) - pom_essentials, goeaobj_vs_orthologs, gene_list_filename_template % ("alb_cer_vs_orth"))
                run_go_on_gene_list((alb_essentials_with_orthologs & pom_essentials) - cer_essentials, goeaobj_vs_orthologs, gene_list_filename_template % ("alb_pom_vs_orth"))
                run_go_on_gene_list((cer_essentials & pom_essentials) - alb_essentials_with_orthologs, goeaobj_vs_orthologs, gene_list_filename_template % ("cer_pom_vs_orth"))
    
    # Write out the group legends:
    write_group_legends(feature_groups, os.path.join(output_folder, "group_legend.csv"))
    
    # Verdict summary:
    with open(os.path.join(output_folder, "essentials_summary.csv"), 'w') as summary_file:
        writer = csv.writer(summary_file)
        writer.writerow(["Group", "Classifier", "RD filter", "PreEvo Hits", "PostEvo Hits", "AUC", "FPR", "TPR", "Threshold", "Essentials", "Sick", "Weirdos", "Non-essentials"])
        for read_depth_filter, dataset in post_analyzed_datasets.iteritems():
            for (classifier_type, group_name) in cls_keys:
                threshold_ix = threshold_ixs[(classifier_type, group_name, read_depth_filter)]
                fpr, tpr, ts = classifiers[(classifier_type, group_name, read_depth_filter)][1]
                cls_auc = auc(fpr, tpr)
                key = "%s-%s-verdict" % (classifier_type, group_name)
                verdicts = [r[key] for r in dataset]
                writer.writerow([group_name, classifier_type, read_depth_filter,
                                 len([h for h in pre_prediction_hits if h["hit_count"] >= read_depth_filter]),
                                 len([h for h in post_prediction_hits if h["hit_count"] >= read_depth_filter]),
                                 cls_auc, fpr[threshold_ix], tpr[threshold_ix], ts[threshold_ix],
                                 verdicts.count("Essential"), verdicts.count("Sick"),
                                 verdicts.count("Weirdo"), verdicts.count("Non-essential")])
    
    # Comparison between all pairs of classifiers:
    stats_keys = ["Essential", "Non-essential", "Weirdo", "Sick"]
    for key1, key2 in sorted(itertools.combinations(cls_keys, 2)):
        (cls_type1, group1) = key1
        (cls_type2, group2) = key2
        verdict_key1 = "%s-%s-verdict" % key1
        verdict_key2 = "%s-%s-verdict" % key2
        for read_depth_filter, dataset in post_analyzed_datasets.iteritems():
            stats = {sk1: {sk2: 0 for sk2 in stats_keys} for sk1 in stats_keys} # Breakdown of key2 by key1
            for record in dataset:
                stats[record[verdict_key1]][record[verdict_key2]] += 1
            with open(os.path.join(output_folder, "comparison_%s-%s_vs_%s-%s.filter-%d.csv" % (cls_type1, group1, cls_type2, group2, read_depth_filter)), 'w') as abs_out:
                writer = csv.writer(abs_out)
                writer.writerow([""] + stats_keys + ["Total", "Agreement"])
                for sk1 in stats_keys:
                    total = sum(stats[sk1].values())
                    writer.writerow([sk1] + [stats[sk1][sk2] for sk2 in stats_keys] + [total, "%.3f" % (stats[sk1][sk1] / float(total))])
                writer.writerow(["Total"] + [sum(stats[sk1][sk2] for sk1 in stats_keys) for sk2 in stats_keys])
                writer.writerow(["Agreement"] + [("%.3f" % (float(stats[sk2][sk2]) / sum(stats[sk1][sk2] for sk1 in stats_keys))) for sk2 in stats_keys])

def draw_venn(title, output_file_path, data, labels):
    plt.figure(figsize=(8,8))
    
    if len(data) == 2:
        venn = venn2
    else:
        venn = venn3
    venn(data, labels)
    plt.title(title)
    plt.savefig(output_file_path,
                transparent=True,
                dpi=300)
    
    plt.close()
    
def write_group_legends(feature_groups, output_file):
    with open(output_file, 'w') as out_file:
        writer = csv.writer(out_file, delimiter=',')
        writer.writerow(["Group name", "Features"])
        for group_name, columns_for_classification in feature_groups:
            writer.writerow([group_name, ", ".join(columns_for_classification)])

def only_orfs(records):
    return (r for r in records if "ORF" in r["feature"].type)

def no_false_negatives(records):
    return (r for r in records if not (r["hits"] == 0 and r["feature"].standard_name not in CER_ESSENTIALS))

def only_consensus_records(records):
    consensus_features = CER_ESSENTIALS | CER_NON_ESSENTIALS
    return (r for r in records if r["feature"].standard_name in consensus_features)

def no_dubious_features(records):
    return (r for r in records if r["feature"].feature_qualifier != "Dubious")

def no_depleted_hit_features(records):
    # If there were no insertions in the feature AND its neighborhood,
    # we can't trust this feature (maybe it was deleted):
    return (r for r in records if r["neighborhood_hits"] > 0 or r["hits"] > 0)

def filter_cer_training_data(cer_records):
    # NB: no_false_negatives - skews the classifiers when training on filtered data
    return no_depleted_hit_features(no_dubious_features(only_orfs(only_consensus_records(no_false_negatives(cer_records)))))   

def _get_sorted_hit_files(folder):
    return sorted(glob.glob(os.path.join(folder, "*_Hits.txt")), key=lambda path: int(re.findall("\d+", path)[-1]))

def _concat_lists(lists):
    result = []
    for l in lists:
        result.extend(l)
    return result

def _take_ixes(seq, ixes):
    return [item for ix, item in enumerate(seq) if ix in ixes]


if __name__ == "__main__":
    output_folder = os.path.join(Shared.get_script_dir(), "output")
    
    all_track_files = glob.glob(os.path.join(Shared.get_script_dir(), "dependencies", "Kornmann", "*.bed"))
    all_track_filenames = [os.path.split(file_path)[-1][:-4] for file_path in all_track_files]
    
    # We cache the hits because the hit reading process involves finding the
    # feature which got hit, for every hit, and that makes reading an O(n log n)
    # operation, which is a little slow. O(n) is better here.
    hit_cache = Shared.get_dependency("Kornmann/cached_sc_track_hits.dat")
    if not os.path.exists(hit_cache):
        all_tracks = [get_hits_from_bed(fname) for fname in all_track_files]
        with open(hit_cache, 'wb') as pickle_file:
            cPickle.dump(all_tracks, pickle_file)
    else:
        with open(hit_cache, 'rb') as pickle_file:
            all_tracks = cPickle.load(pickle_file)
    
    wt2_hits = all_tracks[all_track_filenames.index("Wild Type 2")]
    
    # Test training from cache:
    wt2_analyzed = SummaryTable.analyze_hits(wt2_hits, cer_db, 10000).values()
    wild_type_2_base_records = list( filter_cer_training_data(wt2_analyzed) )
    test_training_data(wild_type_2_base_records,
                       CER_ESSENTIALS,
                       os.path.join(output_folder, "predictions/training/Wild Type 2"))

    pre_hits = SummaryTable.read_hit_files(_get_sorted_hit_files(Shared.get_dependency("albicans/experiment data/pre evo/q20m2")))
    post_hits = SummaryTable.read_hit_files(_get_sorted_hit_files(Shared.get_dependency("albicans/experiment data/post evo/q20m2")))


    classify_albicans(_concat_lists(_take_ixes(pre_hits, (3-1, 7-1, 11-1))),
                      _concat_lists(_take_ixes(post_hits, (3-1, 7-1, 11-1))),
                      wt2_hits,
                      CER_ESSENTIALS,
                      os.path.join(output_folder, "predictions/high cov - wt2 - fpr 0.075"),
                      target_fpr=0.075)
     
    classify_albicans(_concat_lists(_take_ixes(pre_hits, (3-1, 7-1, 11-1))),
                      _concat_lists(_take_ixes(post_hits, (3-1, 7-1, 11-1))),
                      wt2_hits,
                      CER_ESSENTIALS,
                      os.path.join(output_folder, "predictions/high cov - wt2 - fpr 0.05"),
                      target_fpr=0.05)