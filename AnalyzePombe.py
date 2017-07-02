import argparse
import csv
import os
import SummaryTable
import GenomicFeatures
import pandas as pd
import Shared
import Classifier
import glob
import cPickle

# TODO: this scipt is currently very quick-and-dirty in implementation - a lot
# of code duplication, etc. Eventually should be refactored.

def read_hit_file(filename, read_depth_filter=1):
    result = []
    
    with open(filename, "r") as in_file:
        reader = csv.reader(in_file)
        reader.next()
        for line in reader:
            chrom, strand, hit_pos, hit_count, gene_name = line
            
            hit_pos = int(hit_pos)
            hit_count = int(hit_count)
            if hit_count < read_depth_filter:
                continue
            
            obj = {"chrom": chrom,
                   "source": strand,
                   "hit_pos": hit_pos,
                   "hit_count": hit_count,
                   "gene_name": gene_name}
            
            result.append(obj)
                
    return result

def enrich_with_pombe(records):
    pom_db = GenomicFeatures.default_pom_db()
    
    viability_table = pd.read_csv(Shared.get_dependency("pombe/FYPOviability.tsv"),
                                  header=None,
                                  delimiter='\t',
                                  names=["pombe standard name", "essentiality"])
    
    ortholog_table = pd.read_csv(Shared.get_dependency("albicans/C_albicans_SC5314_S_pombe_orthologs.txt"),
                                 skiprows=8,
                                 delimiter='\t',
                                 header=None,
                                 usecols=['albicans standard name', 'pombe standard name'],
                                 names=['albicans standard name', 'albicans common name', 'albicans alb_db id',
                                        'pombe standard name', 'pombe common name', 'pombe alb_db id'])
    
    # TODO: we probably don't want to use the hit table, though the InParanoid
    # table is very stringent.
    best_hit_table = pd.read_csv(Shared.get_dependency("albicans/C_albicans_SC5314_S_pombe_best_hits.txt"),
                                 skiprows=8,
                                 delimiter='\t',
                                 header=None,
                                 usecols=['albicans standard name', 'pombe standard name'],
                                 names=['albicans standard name', 'albicans common name', 'albicans alb_db id',
                                        'pombe standard name', 'pombe common name', 'pombe alb_db id'])
    
    ortholog_table = pd.concat([ortholog_table, best_hit_table])
     
    joined_table = pd.merge(ortholog_table, viability_table, on="pombe standard name")
     
    essentiality_map = {"viable": "No", "inviable": "Yes", "condition-dependent": "Conditional"}
    for record in records:
        ortholog_row = joined_table[joined_table["pombe standard name"] == record["feature"].standard_name]
        if ortholog_row.empty:
            ortholog_name = ""
            ortholog_essentiality = ""
        else:
            ortholog_name = pom_db.get_feature_by_name(ortholog_row["pombe standard name"].iloc[0]).name 
            ortholog_essentiality = essentiality_map.get(ortholog_row["essentiality"].iloc[0], "?")
         
        record["pombe_ortholog"] = ortholog_name
        record["essential_in_pombe"] = ortholog_essentiality

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--hits-file", required=True,
                        help="The hit file to analyze.")
    parser.add_argument("--rdf", default=1, type=int,
                        help="Minimum read depth filter per insertion to include in the analysis.")
    
    args = parser.parse_args()
    hits_file = args.hits_file
    rdf = args.rdf
    hits_filename = os.path.splitext(os.path.split(hits_file)[1])[0]
    
    out_file = os.path.splitext(hits_file)[0] + ".rdf_%d.analyzed.csv" % rdf
    hits = read_hit_file(hits_file, rdf)
    
    pom_db = GenomicFeatures.default_pom_db()
    
    analyzed_hits = SummaryTable.analyze_hits(hits, pom_db)
    records = analyzed_hits.values()
    enrich_with_pombe(records)
    
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
            "field_name": "essential_in_pombe",
            "csv_name": "Essential in Sp",
        },
        
        {
            "field_name": "RF-G3",
            "csv_name": "RF-G3",
        },
        
        {
            "field_name": "hits",
            "csv_name": "Hits",
            "format": "%d"
        },
        
        {
            "field_name": "reads",
            "csv_name": "Reads",
            "format": "%d"
        },
        
        {
            "field_name": "feature",
            "csv_name": "Length",
            "format": lambda f: len(f)
        },
        
        {
            "field_name": "neighborhood_index",
            "csv_name": "Neighborhood index",
            "format": "%.3f"
        },
        
        {
            "field_name": "upstream_hits_100",
            "csv_name": "100 bp upstream hits"
        },
        
        {
            "field_name": "max_free_region",
            "csv_name": "Max free region"
        },
                   
        {
            "field_name": "freedom_index",
            "csv_name": "Freedom index",
            "format": "%.3f"
        },
        
        {
            "field_name": "kornmann_domain_index",
            "csv_name": "Kornmann DI",
            "format": "%.3f"
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
    
    output_folder = "."
    
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
    cer_db = GenomicFeatures.default_cer_db()
    wt2_analyzed = SummaryTable.analyze_hits(wt2_hits, cer_db, 10000).values()
    wild_type_2_base_records = list( Classifier.filter_cer_training_data(wt2_analyzed) )

    CER_ESSENTIALS, CER_NON_ESSENTIALS = SummaryTable.get_cerevisiae_essentials()
    target_fpr=0.075
    
    from sklearn.metrics import roc_curve, auc
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import RandomForestClassifier

    classifier_factory = {
        "LR": lambda: LogisticRegression(),
        "RF": lambda: RandomForestClassifier(n_estimators=100, random_state=0)
    }
    
    feature_groups = (
        ("G3", ("neighborhood_index", "length", "hits", "freedom_index", "upstream_hits_100")),
    )
    
    classifiers = Classifier.train_classifiers(
        wild_type_2_base_records,
        CER_ESSENTIALS,
        feature_groups,
        classifier_factory
    )
        
    for r in records:
        r["RF-G3"] = "N/A"
    
    print "There are %d false non-essentials" % len([r for r in records if r["essential_in_pombe"] == "No" and r["hits"] == 0])
    known_pombe_genes = [record for record in records
                         if record["essential_in_pombe"] in ("Yes", "No")
#                          and not (record["essential_in_pombe"] == "No" and record["hits"] == 0)
                         ]
    annotations = [1 if record["essential_in_pombe"] == "Yes" else 0
                   for record in known_pombe_genes]
    
    # TODO: quite ugly, fix this:
    classifiers = {(k[0], k[1], None): v for k, v in classifiers.items()}
    
    Classifier.predict_records(
        known_pombe_genes,
        feature_groups,
        classifiers
    )
    
    SummaryTable.write_data_to_csv(records, cols_config, out_file)
    
    classifier_name = "RF"
    fpr, tpr, thresholds = roc_curve(annotations, [record["RF-G3"] for record in known_pombe_genes])
    classifier_auc = auc(fpr, tpr)
    
    with open(os.path.join(output_folder, hits_filename + ".fpr_table.rdf_%d.csv" % rdf), "w") as out_file:
        writer = csv.writer(out_file)
        writer.writerow(["FPR", "TPR", "Threshold"])
        for a, b, c in zip(fpr, tpr, thresholds):
            writer.writerow([a, b, c])
    
    import matplotlib.pyplot as plt
    
    plt.plot(fpr, tpr, label="%s (AUC:%.2f)" % (classifier_name, auc(fpr, tpr)))
    plt.legend(loc="lower right")
    plt.savefig(os.path.join(output_folder, '%s.ROC.rdf_%d.pdf' % (hits_filename, rdf)),
                transparent=True,
                dpi=300)
    plt.close()
    