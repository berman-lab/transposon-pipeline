import os
import csv
import GeneFeatures
from GeneFeatures import cer_db, alb_db
import argparse
import glob
import scipy.stats
import numpy as np
import math
import pandas as pd
import Shared

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle

def read_hit_files(files, read_depth_filter=1):
    """Read in the list of hits files.
    
    Attributes
    ----------
    read_depth_filter : int
        The read depth below which insertion events will be ignored.
    """

    return [read_hit_file(f, read_depth_filter) for f in files]

def read_hit_file(filename, read_depth_filter=1):
    """Read in the given hit file.
    
    Attributes
    ----------
    read_depth_filter : int
        The read depth below which insertion events will be ignored.
    
    Returns
    -------
    list of dicts of str to object
        Every hit is converted to a dict, specifying its properties according
        to the original table columns. Numerical values are converted as
        necessary.
    """
    
    result = []
    
    with open(filename, "r") as in_file:
            in_file.next()
            for line in in_file:
                chrom, source, up_feature_type, up_feature_name, up_gene_name, \
                       up_feature_dist, down_feature_type, down_feature_name, \
                       down_gene_name, down_feature_dist, ig_type, \
                       hit_pos, hit_count = line.split("\t")
                
                hit_count = int(hit_count)
                if hit_count < read_depth_filter:
                    continue       
                
                if "ORF" in ig_type:
                    gene_name = up_gene_name if up_gene_name != "nan" \
                            else up_feature_name
                    assert ig_type in ("ORF(W)-ORF(W)", "ORF(C)-ORF(C)") and \
                           up_feature_name == down_feature_name and \
                           up_gene_name == down_gene_name and \
                           up_feature_dist == down_feature_dist == "0"
                else:
                    gene_name = "nan"
                
                obj = {"chrom": chrom,
                       "source": source,
                       # NB: these fields are not used, but take up a lot of
                       # space. Ignore them for now.
#                        "up_feature_type": up_feature_type,
#                        "up_feature_name": up_feature_name,
#                        "up_gene_name": up_gene_name,
#                        "up_feature_dist": up_feature_dist,
#                        "down_feature_type": down_feature_type,
#                        "down_feature_name": down_feature_name,
#                        "down_gene_name": down_gene_name,
#                        "down_feature_dist": down_feature_dist,
#                        "ig_type": ig_type,
                       "hit_pos": int(hit_pos),
                       "hit_count": hit_count,
                       "gene_name": gene_name}
                
                result.append(obj)
                
    return result

TOTAL_HITS = "Total Hits"
TOTAL_READS = "Total Reads"
AVG_READS_PER_HIT = "Mean Reads Per Hit"
ORF_HITS = "Hits in ORFs"
ANNOTATED_FEATURE_HITS = "Hits in Genomic Features"
PER_ANNOTATED_FEATURE_HITS = "% of hits in features"
INTERGENIC_HITS = "Intergenic Hits"
PER_INTERGENIC_HITS = "% of intergenic hits"
FEATURES_HIT = "No. of Features Hit"
PER_FEATURES_HIT = "% of features hit"
AVG_HITS_IN_FEATURE = "Mean hits per feature"
AVG_READS_IN_FEATURE = "Mean reads per feature"
AVG_READS_IN_FEATURE_HIT = "Mean reads per hit in feature"
READS_IN_FEATURES = "Total reads in features"

ALL_STATS = [TOTAL_READS, TOTAL_HITS, PER_ANNOTATED_FEATURE_HITS, PER_INTERGENIC_HITS,
             PER_FEATURES_HIT, AVG_HITS_IN_FEATURE, AVG_READS_IN_FEATURE, 
             AVG_READS_PER_HIT, AVG_READS_IN_FEATURE_HIT]

def get_statistics(dataset, feature_db):
    """Get general statistics about hit files.
    
    Returns
    -------
    dict
        A dict containing the computed stats.
    """
    
    result = {}
    
    result[TOTAL_HITS] = len(dataset)
    result[TOTAL_READS] = sum(obj["hit_count"] for obj in dataset)
    result[AVG_READS_PER_HIT] = result[TOTAL_READS] / result[TOTAL_HITS]
    result[ORF_HITS] = 0
    result[ANNOTATED_FEATURE_HITS] = 0
    result[INTERGENIC_HITS] = 0
    result[READS_IN_FEATURES] = 0
    
    features_hit = set()
    for hit in dataset:
        features = feature_db.get_features_at_location(hit["chrom"], hit["hit_pos"])
        if len(features) == 0:
            continue
        result[ANNOTATED_FEATURE_HITS] += 1
        result[READS_IN_FEATURES] += hit["hit_count"]
        features_hit.update(set( f.standard_name for f in features ))
        for f in features:
            if f.is_orf:
                result[ORF_HITS] += 1
                break
    
    result[INTERGENIC_HITS] = result[TOTAL_HITS] - result[ANNOTATED_FEATURE_HITS]
    result[FEATURES_HIT] = len(features_hit)
    result[PER_INTERGENIC_HITS] = "%.2f%%" % (result[INTERGENIC_HITS] * 100.0 /  result[TOTAL_HITS])
    result[PER_ANNOTATED_FEATURE_HITS] = "%.2f%%" % (result[ANNOTATED_FEATURE_HITS] * 100.0 /  result[TOTAL_HITS])
    result[PER_FEATURES_HIT] = "%.2f%%" % (result[FEATURES_HIT] * 100.0 / len(alb_db.get_all_features()))
    result[AVG_HITS_IN_FEATURE] = "%.1f" % (result[ANNOTATED_FEATURE_HITS] * 1.0 / result[FEATURES_HIT])
    result[AVG_READS_IN_FEATURE] = result[READS_IN_FEATURES] / result[FEATURES_HIT]
    result[AVG_READS_IN_FEATURE_HIT] = result[READS_IN_FEATURES] / result[ANNOTATED_FEATURE_HITS]
    
    return result

def analyze_hits(dataset, feature_db, neighborhood_window_size=10000):
    """Analyze a hit dataset into a feature dataset.
    
    Attributes
    ----------
    dataset : list of dict
        The hit list, requires that each hit have a "hit_pos", "hit_count" and
        "gene_name" field.
    feature_db : `GeneFeatures._FeatureDB`
        The feature database to use.
    neighborhood_window_size : int
        The window size for computing the neighborhood index.
    
    Returns
    -------
    dict of dicts
        A map of standard feature names to its all_analyzed results. The results
        themselves are represented as a dict of str to object. For example:
        
        {
            "C2_":
                {
                    "feature": the_feature_object,
                    "hits": 24,
                    "reads": 5064,
                    ... 
                }
            ...,
        }
    """
    
    log2 = lambda v: math.log(v, 2)
    
    result = {}
    
    total_reads = sum(h["hit_count"] for h in dataset)
    total_reads_log = log2(total_reads)
    
    chroms = set(h["chrom"] for h in dataset)
    for chrom in chroms:
        hits = [h for h in dataset if h["chrom"] == chrom]
        
        chrom_len = max(len(feature_db[chrom]), max(h["hit_pos"] for h in hits))
        # For ease of analysis downstream, we separate the hits that fell within
        # features from those that fell outside.
        hits_in_features = np.zeros((chrom_len + 1,), dtype=np.int)
        hits_outside_features = np.zeros((chrom_len + 1,), dtype=np.int)
        reads_across_chrom = np.zeros((chrom_len + 1,), dtype=np.int)
        
        for hit in hits:
            reads_across_chrom[hit["hit_pos"]] += hit["hit_count"]
            # TODO: why do we cap the hits at 2? What happens with pooled
            # hit all_hits?
            if hit["gene_name"] == "nan":
                if hits_outside_features[hit["hit_pos"]] < 2:
                    hits_outside_features[hit["hit_pos"]] += 1 
            else:
                if hits_in_features[hit["hit_pos"]] < 2:
                    hits_in_features[hit["hit_pos"]] += 1
        
        records = {} # The per-feature records for this chromosome.

        for feature in feature_db[chrom]:
            hits_in_feature = hits_in_features[feature.start:feature.stop+1]
            hits_in_feature_count = hits_in_feature.sum()
            reads_in_feature = reads_across_chrom[feature.start:feature.stop+1].sum()
            
            # The borders of the neighborhood window:
            window_start = max(1, feature.start - neighborhood_window_size)
            window_end = min(chrom_len, feature.stop + neighborhood_window_size)
            
            hits_outside_feature = hits_outside_features[window_start:window_end+1]
            hits_outside_feature_count = hits_outside_feature.sum()
            
            intergenic_region = feature_db.get_interfeature_range(chrom, (window_start, window_end))
            insertion_index = float(hits_in_feature_count) / len(feature)
            
            # There are some scenarios where large regions have no hits whatsoever:
            #   * The dataset is very low coverage.
            #   * A read depth filter was applied, and it greatly decreased coverage.
            #   * A region is simply deleted.
            # Computing the neighborhood index would thus result in a division-by-zero
            # error. To prevent this, we assign the neighborhood index to be 0,
            # however it's up to the downstream tools to filter this. It's suggested
            # to ignore features that had no hits in their neighborhood, as they're
            # probably not informative and will increase error rates.
            if hits_outside_feature_count == 0:
                neighborhood_index = 0
            else:
                neighborhood_index = insertion_index / (float(hits_outside_feature_count) / intergenic_region.coverage)
            
            # TODO: We assume that the 100 bp upstream is always intergenic.
            if feature.strand == 'W':
                upstream_slice_100 = slice(max(1, feature.start - 100), feature.start)
                upstream_slice_50 = slice(max(1, feature.start - 50), feature.start)
            elif feature.strand == 'C':
                upstream_slice_100 = slice(feature.stop + 1, min(chrom_len, feature.stop + 101))
                upstream_slice_50 = slice(feature.stop + 1, min(chrom_len, feature.stop + 51))
            else:
                # TODO: does this happen? Can this happen? What should we do about it?
                # This can happen in cerevisiae ARS regions.
                assert False
                upstream_slice_100 = slice(0, 0)
                upstream_slice_50 = slice(0, 0)
            
            upstream_region_100 = hits_outside_features[upstream_slice_100]
            upstream_region_50 = hits_outside_features[upstream_slice_50]
            
            # Compute the longest area without any hits:
            hit_ixes = [0] + list(np.where(hits_in_feature > 0)[0]+1) + [len(feature) + 1]
            max_free_region = max(right - left for left, right in zip(hit_ixes, hit_ixes[1:])) - 1
            
            upstream_50_hits = upstream_region_50.sum()
            upstream_100_hits = upstream_region_100.sum()
            
            records[feature.standard_name] = {
                "feature": feature,
                "length": len(feature),
                "hits": hits_in_feature_count,
                "reads": reads_in_feature,
                # Can be tested downstream for zero:
                "neighborhood_hits":  hits_outside_feature_count,
                "insertion_index": insertion_index,
                "neighborhood_index": neighborhood_index,
                "upstream_hits_50": upstream_50_hits,
                "upstream_hits_100": upstream_100_hits,
                "max_free_region": max_free_region,
                "freedom_index": float(max_free_region) / len(feature),
                # The value is singular, to get the coefficient we subtract the pre from the post later on:
                "s_value": log2(reads_in_feature + 1) - total_reads_log, # Add 1 so as to not get a log 0
                # Note: the positions are relative to the genome, NOT the gene:
                "hit_locations": [ix+1 for (ix, hit) in enumerate(hits_in_feature) if hit > 0],
            }
            
        result.update(records)
    
    return result

def write_analyzed_alb_records(records, output_file):
    """Write the analyzed albicans records with a default column configuration
    into a csv file.
    
    Parameters
    ----------
    records : list of dicts
        A list of analyzed records.
    output_file : str
        The path to the output file.
    """
    
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
            "field_name": "essential_in_cerevisiae",
            "csv_name": "Essential in Sc",
        },
        
        {
            "field_name": "cer_synthetic_lethal",
            "csv_name": "Sc synthetic lethals",
        },
        
        {
            "field_name": "cer_fitness",
            "csv_name": "Sc fitness",
        },
        
        {
            "field_name": "essential_in_albicans",
            "csv_name": "Essential in albicans",
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
    
    enrich_alb_records(records)
    
    write_data_to_csv(records, cols_config, output_file)

def write_data_to_csv(data, cols_config, output_filename):
    """Write a list of dicts representing data objects into a CSV file using a
    specified configuration.
    
    The configuration allows the standard C-style string formatting as well as
    passing a custom function for data transformation.
    
    Note
    ----
    
    This is a legacy method. Consider using pandas instead.
    
    Parameters
    ----------
    data : iterable of dicts of str to object
        An iterable of data objects.
    cols_config : list of dicts of str to object 
        An example of a configuration:
        [
            {
                field_name: "field_name",    # Mandatory, the field name in the data object dicts.
                format: "%.3f",              # Optional, default is %s. Can be a callable - in which case, the callable will be executed on the field.
                csv_name: "Column name",     # Mandatory, must be unique.
                sort_by: True,               # Optional, default is not to sort by any column. Can only be set once. Cannot be set with a callable for the format.
            },
            ...
        ]
    output_filename : str
        The name of the output file.
    """
    
    field_col = "field_name"
    format_col = "format"
    csv_name_col = "csv_name"
    sort_by_col = "sort_by"
    
    # Filter out field names that don't exist:
    filtered_cols_config = []
    first_record = data[0]
    for col in cols_config:
        if col[field_col] not in first_record:
            print "WARNING: %s not in records to print!" % col[field_col]
            continue
        filtered_cols_config.append(col)
    cols_config = filtered_cols_config
    
    # Find sort column:
    sort_field = None
    for col in cols_config:
        if col.get(sort_by_col, False):
            sort_field = col[field_col]
            break
    
    # Sort by the sort column, if found:
    if sort_field is not None:
        ordered_dataset = sorted(data, key=lambda record: record[sort_field])
    else:
        ordered_dataset = data # Default ordering, really
    
    # Cache the formats:
    format_cache = {}
    for col in cols_config:
        col_name = col[csv_name_col]
        if format_col not in col:
            format_cache[col_name] = lambda s: "%s" % s
        elif callable(col[format_col]):
            format_cache[col_name] = col[format_col]
        else:
            format_cache[col_name] = lambda s: col[format_col] % s
    
    # Write the output CSV:
    with open(output_filename, 'w') as out_file:
        writer = csv.writer(out_file, delimiter=',')
        
        # Write out the header:
        writer.writerow([col[csv_name_col] for col in cols_config])
        
        # Write the entire dataset:
        for record in ordered_dataset:
            writer.writerow([format_cache[col[csv_name_col]](record[col[field_col]]) for col in cols_config])

def get_cerevisiae_essentials():
    """Get consensus essentials and non-essentails in cerevisiae.
    
    Consensus is determined by literature agreement - all annotations are
    required to state one or the other. If disagreement exists, the gene
    isn't included in any dataset.
    
    Returns
    -------
    (set, set)
        A pair sets, denoting the essential and non-essential genes, using
        their standard names.
    """
    
    viable_filepath = Shared.get_dependency("cerevisiae/cerevisiae_viable_annotations.txt")
    inviable_filepath = Shared.get_dependency("cerevisiae/cerevisiae_inviable_annotations.txt")
    
    viable_table = pd.read_csv(viable_filepath, skiprows=8, delimiter="\t")
    inviable_table = pd.read_csv(inviable_filepath, skiprows=8, delimiter="\t")
    
    annotated_as_viable = set(cer_db.get_feature_by_name(f) for f in set(viable_table["Gene"])) - set([None])
    annotated_as_inviable = set(cer_db.get_feature_by_name(f) for f in set(inviable_table["Gene"])) - set([None])
    
    consensus_viable = annotated_as_viable - annotated_as_inviable
    consensus_inviable = annotated_as_inviable - annotated_as_viable
    
    # TODO: the dubious genes shouldn't be filtered here.
    consensus_viable_orfs = [f for f in consensus_viable if f.is_orf and f.feature_qualifier != "Dubious"]
    consensus_inviable_orfs = [f for f in consensus_inviable if f.is_orf and f.feature_qualifier != "Dubious"]
    
    return (set(f.standard_name for f in consensus_inviable_orfs),
            set(f.standard_name for f in consensus_viable_orfs))

def get_alb_coverage():
    """Get a measure of uniqueness of each albicans gene.
    
    Unique coverage is defined by the relative area of the gene that can be
    covered by 100-150 bp reads with a mapping quality of >= 20. This was
    derived from an NGS analysis of the parental haploid.
    
    Returns
    -------
    dict of str to float
        A map of standard feature names to a 0-1 float denoting the relative
        uniqueness of coverage. 
    """
    
    unique_coverage = {}
    with open(Shared.get_dependency("albicans/unique_coverage.csv"), 'r') as in_file:
        reader = csv.reader(in_file, delimiter=',')
        for line in reader:
            feature_name, percent_covered = line
            unique_coverage[feature_name] = float(percent_covered) / 100
        
    return unique_coverage

def enrich_alb_records(records):
    """Enrich albicans records with 3-rd party data.
    
    3-rd party data include orthologs, essentiality from literature, etc."""
    
    # Pombe:
    enrich_with_pombe(records)

    # Cerevisiae:
    cer_essentials, cer_non_essentials = get_cerevisiae_essentials()

    # TODO: document the origin of the synthetic lethal dataset.
    cer_synthetic_lethals = {}
    with open(Shared.get_dependency("cerevisiae/duplicatesSl_011116.txt"), 'r') as in_file:
        in_file.readline()
        for line in in_file:
            f1_name, f2_name, _, is_double_lethal = line.split()
            
            f1 = cer_db.get_feature_by_name(f1_name)
            f2 = cer_db.get_feature_by_name(f2_name)
            
            f1_std_name = f1.standard_name if f1 is not None else ""
            f2_std_name = f2.standard_name if f2 is not None else ""
            
            cer_synthetic_lethals[f1_std_name] = cer_synthetic_lethals[f2_std_name] = \
                {"1": "Yes", "0": "No"}.get(is_double_lethal, is_double_lethal)
    
    cer_fitness = {}
    # TODO: document the origin of the fitness dataset.
    with open(Shared.get_dependency("cerevisiae/neFitnessStandard.txt"), 'r') as in_file:
        for line in in_file:
            fname, fitness = line.split()
            feature = cer_db.get_feature_by_name(fname)
            if feature is None:
                continue 
            cer_fitness[cer_db.get_feature_by_name(fname).standard_name] = fitness
            
    # Albicans:
    
    # All liteartue phenotypes, no pre-filtering:
    
    alb_phenotype_table_path = Shared.get_dependency("albicans/C_albicans_phenotype_data.tab")
    alb_phenotype_table = pd.read_csv(
        alb_phenotype_table_path,
        delimiter='\t',
        header=None,
        usecols=['Feature Name', 'Phenotype'],
        names=['Feature Name', 'Feature Type', 'Gene Name',
               'CGDID', 'Reference', 'Experiment Type', 'Mutant Type',
               'Allele', 'Strain background', 'Phenotype', 'Chemical', 'Condition',
               'Details', 'Reporter', 'Anatomical Structure', 'Virulence Model', 'Species']
    )
    
    alb_pheno_essentials = set(alb_phenotype_table[alb_phenotype_table["Phenotype"] == "inviable"]["Feature Name"])
    alb_pheno_non_essentials  = set(alb_phenotype_table[alb_phenotype_table["Phenotype"] == "viable"]["Feature Name"])
    
    alb_db = GeneFeatures.alb_db
    alb_pheno_essential_features = set(f.standard_name for f in filter(None, (alb_db.get_feature_by_name(f) for f in alb_pheno_essentials)))
    alb_pheno_non_essential_features = set(f.standard_name for f in  filter(None, (alb_db.get_feature_by_name(f) for f in alb_pheno_non_essentials)))
    
    alb_pheno_toss_up = alb_pheno_essential_features & alb_pheno_non_essential_features
    alb_pheno_essential_consensus = alb_pheno_essential_features - alb_pheno_non_essential_features
    alb_pheno_non_essential_consensus = alb_pheno_non_essential_features - alb_pheno_essential_features
    
    grace_table_path = Shared.get_dependency("albicans/ncomms7741-s2.xls")
    grace_table = pd.read_excel(grace_table_path, sheetname=1, #"Essentiality scores",
                                skiprows=14, header=None,
                                names=["orf19 name", "Common", "Description", "Plate", "Position",
                                       "tet growth", "5-FOA excision", "dox growth", "Essentiality Concordance  (Y/N)",
                                       "Essential/Not essential", "S. cerevisiae homolog", "S. cerevisiae KO phenotype"])
    
    # GRACE:
    
    roemer_grace_essentials = set()
    roemer_grace_non_essentials = set()
    omeara_grace_essentials = set()
    omeara_grace_non_essentials = set()
    
    # Couldn't figure out how to make pandas constrain types in a column, so we
    # have to handle all sorts of weirdness, including nan, and make it as
    # robust as possible.
    def _extract_tet_growth(tet_growth):
        str_tet_growth = str(tet_growth)
        if "+" not in str_tet_growth:
            return None
        return float(str_tet_growth.replace("+", ""))
    
    def _extract_foa_excision(foa_excision):
        if foa_excision not in ("Yes", "No"):
            return None
        return foa_excision
        
    for _ix, line in grace_table.iterrows():
        orf19name = line["orf19 name"]
        feature = alb_db.get_feature_by_name(orf19name)
        if feature is None:
            continue
        
        omeara_essential = line["Essential/Not essential"]
        if omeara_essential == 'E':
            omeara_grace_essentials.add(feature.standard_name)
        elif omeara_essential == 'N':
            omeara_grace_non_essentials.add(feature.standard_name) 
        
        tet_growth = _extract_tet_growth(line["tet growth"])
        foa_excision = _extract_foa_excision(line["5-FOA excision"])
        
        
        if tet_growth is None and foa_excision is None:
            continue
        
        if (tet_growth >= 3) or foa_excision == "Yes":
            roemer_grace_essentials.add(feature.standard_name)
        else:
            roemer_grace_non_essentials.add(feature.standard_name)
    
    unique_coverage = get_alb_coverage()
    
    for record in records:
        feature = record["feature"]
        if feature.cerevisiae_orthologs:
            # TODO: according to albicans feature file, a gene may have more than one ortholog, but in practice this doesn't happen.
            assert len(feature.cerevisiae_orthologs) == 1
            cer_ortholog = cer_db.get_feature_by_name(iter(feature.cerevisiae_orthologs).next()).standard_name
            
            record["essential_in_cerevisiae"] = "Yes" if cer_ortholog in cer_essentials else \
                                                "No" if cer_ortholog in cer_non_essentials else "?"
             
            record["cer_synthetic_lethal"] = cer_synthetic_lethals.get(cer_ortholog, "")
            record["cer_fitness"] = cer_fitness.get(cer_ortholog, "")
        else:
            record["essential_in_cerevisiae"] = ""
            record["cer_synthetic_lethal"] = ""
            record["cer_fitness"] = ""
        record["unique_coverage"] = unique_coverage.get(feature.standard_name, 1)
        record["essential_in_albicans"] = "Yes" if feature.standard_name in alb_pheno_essential_consensus else \
            "No" if feature.standard_name in alb_pheno_non_essential_consensus else \
            "?" if feature.standard_name in alb_pheno_toss_up else ""
        record["essential_in_albicans_grace_roemer"] = "Yes" if feature.standard_name in roemer_grace_essentials else \
            "No" if feature.standard_name in roemer_grace_non_essentials else ""
        record["essential_in_albicans_grace_omeara"] = "Yes" if feature.standard_name in omeara_grace_essentials else \
            "No" if feature.standard_name in omeara_grace_non_essentials else ""

def enrich_with_pombe(records):
    """Add pombe orthologs and essentiality annotations to albicans records."""
    
    pombe_genes = GeneFeatures.get_pombe_genes()
    
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
        ortholog_row = joined_table[joined_table["albicans standard name"] == record["feature"].standard_name]
        if ortholog_row.empty:
            ortholog_name = ""
            ortholog_essentiality = ""
        else:
            ortholog_name = pombe_genes.get(ortholog_row["pombe standard name"].iloc[0], ortholog_row["pombe standard name"].iloc[0]) 
            ortholog_essentiality = essentiality_map.get(ortholog_row["essentiality"].iloc[0], "?")
         
        record["pombe_ortholog"] = ortholog_name
        record["essential_in_pombe"] = ortholog_essentiality
        
def compare_insertion_and_neighborhood((analysis_labels, analyses), output_file):
    """Compare the insertion and the neighborhood indices with Pearson and
    Spearman correlations.
    
    Parameters
    ----------
    analysis_labels : sequence
        The labels for each analysis dataset.
    analyses : sequence
        A sequence of analyzed datasets (as returned by `analyze_hits`).
    output_file : str
        The path to the output file.
    """
    
    with open(output_file, 'w') as out_file:
        for label, analysis in zip(analysis_labels, analyses):
            s1 = np.array([r["insertion_index"] for r in analysis.values()])
            s2 = np.array([r["neighborhood_index"] for r in analysis.values()])
            pearson = scipy.stats.pearsonr(s1, s2)
            spearman = scipy.stats.spearmanr(s1, s2)
            out_file.write("%s\nPearson: %s\nSpearman: %s\n" % (label, pearson, spearman))    

def draw_histogram_of_analysis(dataset, output_file_prefix):
    # These parameters should perhaps be refactored to be part of the function
    # signature. 
    for field_name, start, stop, bins, ylim in (("neighborhood_index", 0, 1.5, 300, 100),):
        # The relative_values list is the source for the histogram. For every
        # feature it will contain a normalized value which will later be
        # converted into a histogram. The normalization is done within the
        # arbitrary [start, stop] range and the number of bins.
        
        # We filter out the outliers with 0 hits because there is a lot of
        # them and they skew the relative_values.
        relative_values = [min(bins, int(bins*(r[field_name] - start)/(stop - start))) for r in dataset.values() if r["hits"] > 0]
        zero_outliers = len([r for r in dataset.values() if r["hits"] == 0])
        right_tail_outliers = len([v for v in relative_values if v == bins])
        with open(output_file_prefix + "_outlier_stats.txt", 'w') as out_file:
            out_file.write(" ".join(str(o) for o in (field_name, zero_outliers, right_tail_outliers)))
        relative_values = [v for v in relative_values if v < bins] # Clip the long-tail end
        # Should we include the 0 outliers?
#         relative_values = [v+1 for v in relative_values] + [0]*zero_outliers; bins += 1
        relative_values += [0]*zero_outliers
        draw_histogram(relative_values,
                       "",
                       "%s.%s.%.3f_%d.png" % (output_file_prefix, field_name, stop, bins),
                       bins=bins,
                       ylim=ylim)
            
#             with open(os.path.join(output_dir, "%s.for_matlab.%.3f_%d.txt" % (fname, stop, bins)), 'w') as out_file:
#                 new_hist = [0 if h is 0 else h-1 for h in relative_values]
#                 for b in range(bins):
#                     out_file.write("%d\n" % len([h for h in new_hist if h == b]))

def draw_histogram(histogram, sample_name, file_name, bins, ylim=300):
    plt.hist(histogram, bins=bins)
    plt.title(sample_name)
    plt.xlabel("Hit bins")
    plt.ylabel("# of features")
#     plt.axes.set_ylim(top=300)
    axes = plt.gca()
    axes.set_ylim([0, ylim])
    axes.set_xlim([0, bins])
     
    plt.savefig(file_name)
     
    plt.close()
    
def draw_chrom_map(data, bin_size, title, y_max, outfile):
    chroms = list(alb_db)
    chroms.sort(key=lambda c: c.name)
    max_chrom_len = max(len(c) for c in chroms)
    
    plt.figure(figsize=(30, 4*len(chroms)))
    plt.gcf().suptitle(title, fontsize=22, fontweight='bold', y=0.92, x=0.2)
    
    ade2 = alb_db.get_feature_by_name("ADE2")
    
    for chrom_ix, chrom in enumerate(chroms):
        chrom_display_name = chrom.name[4:8].capitalize()
        
        ax = plt.subplot2grid((len(chroms), 100),
                              (chrom_ix, 0),
                              colspan=len(chrom)*100/max_chrom_len)
        
        #adding the last position as 0/ 1 so the plot won't be shorter    
        ax.bar(range(len(data[chrom.name])), data[chrom.name])
        ax.set_ylabel(chrom_display_name, fontsize = 20, fontweight = 'bold')
        ax.set_ylim([0, y_max])
        ax.set_xlim([0, len(data[chrom.name])])
        ax.tick_params(axis='both', which='major', labelsize=14)
    
        patches = {
            "centromere": "red",
            "repeat_region": "green",
            "long_terminal_repeat": "brown",
            "tRNA|Verified": "yellow",
            "tRNA|Uncharacterized": "yellow",
            "retrotransposon": "black"
        }
        for feature in chrom.get_features():
            if feature.type not in patches:
                continue
            ax.add_patch(Rectangle((float(feature.start)/bin_size, 0),
                                   width=float(len(feature))/bin_size,
                                   height=-y_max/10.0,
                                   fc=patches[feature.type],
                                   ec=patches[feature.type],
                                   clip_on=False))
            
        if chrom.name == ade2.chromosome:
            ax.add_patch(Rectangle((ade2.start/bin_size, 0),
                                   width=1,
                                   height=-y_max/10.0,
                                   fc="purple",
                                   ec="purple",
                                   clip_on=False))
            
    plt.savefig(outfile)
    plt.close()

def write_hits_into_bed(target_file, hits):
    """Write a given hit list into a BED-format file.
    
    This is useful for visualizing the hits in genome browsers that support
    the BED format, such as GBrowse."""
    
    with open(target_file, 'w') as bed_file:
        bed_file.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\n")
        for hit in sorted(hits, key=lambda o: (o["chrom"], o["hit_pos"])):
            # NB: BED locations are zero-based, while ours are one-based.
            bed_file.write("%s\t%d\t%d\t.\t%d\t%s\n" %
                           (hit["chrom"], hit["hit_pos"] - 1, hit["hit_pos"], hit["hit_count"],
                            {'W': '+', 'C': '-'}[hit["source"]]))

def write_analyzed_hits_into_bed_proteome(target_file, records):
    """Write a given hit list into a BED-format file that represents protein.
    
    Note
    ----
    
    This function ignores introns - genes with introns will not be displayed correctly!
    """
    
    with open(target_file, 'w') as bed_file:
        bed_file.write("#orf\torfStart\torfEnd\n")
        for record in sorted(records, key=lambda r: r["feature"].standard_name):
            feature = record["feature"]
            aa_hits = sorted(set( ((h-1)/3)+1 for h in record["hit_locations"] ))
            if feature.strand == 'C':
                # As the hit locations are relative to the genome, that is the
                # Watson strand, we have to reverse them if the feature is
                # translated from the Crick strand. 
                feature_len = len(feature) / 3
                aa_hits = [feature_len - h + 1 for h in aa_hits]
            for aa_hit in aa_hits:
                bed_file.write("%s\t%d\t%d\n" % (record["feature"].standard_name, aa_hit, aa_hit+1))

def make_hit_map(hits, bin_size):
    result = {chrom.name: [0]*(len(chrom)/bin_size + 1) for chrom in alb_db}
    
    for hit in hits:
        result[hit["chrom"]][hit["hit_pos"]/bin_size] += 1
    
    return result

def make_read_map(hits, bin_size):
    result = {chrom.name: [0]*(len(chrom)/bin_size + 1) for chrom in alb_db}
    
    for hit in hits:
        result[hit["chrom"]][hit["hit_pos"]/bin_size] += hit["hit_count"]
    
    return result

def transform_read_map_to_log(read_map, base):
    for data in read_map.values():
        for ix in range(len(data)):
            data[ix] = math.log(data[ix], base) if data[ix] > 0 else 0
    
    return read_map

if __name__ == "__main__":    
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-dir", required=True)
    parser.add_argument("-o", "--output-dir", required=True)
    parser.add_argument("-f", "--read-depth-filter", type=int, default=1)
    args = parser.parse_args()
    
    input_dir = args.input_dir
    output_dir = args.output_dir
    read_depth_filter = args.read_depth_filter
    
    Shared.make_dir(output_dir)
    
    input_file_paths = glob.glob(os.path.join(input_dir, "*_Hits.txt"))
    input_filenames = [os.path.split(file_path)[-1][:-9] for file_path in input_file_paths]
    
    all_hits = read_hit_files(input_file_paths, read_depth_filter)
    
    all_analyzed = [analyze_hits(dataset, alb_db, 10000) for dataset in all_hits]
    
    bin_size = 10000
    for fname, hits in zip(input_filenames, all_hits):
        hit_map = make_hit_map(hits, bin_size)
        read_map = make_read_map(hits, bin_size)
        # NB: using a different log base doesn't affect the final figure much,
        # as the relationships between the values stay the same no matter what
        # base is used. 
        log_read_map = transform_read_map_to_log(make_read_map(hits, bin_size), 10)
        draw_chrom_map(hit_map, bin_size, "Hits/10000 bps", 500, os.path.join(output_dir, "hit_map.%s.png" % fname))
        draw_chrom_map(log_read_map, bin_size, "Log10 reads/10000 bps", 7, os.path.join(output_dir, "log10_read_map.%s.png" % fname))
        draw_chrom_map(read_map, bin_size, "Reads/10000 bps", 100000, os.path.join(output_dir, "read_map.%s.png" % fname))
    
    # Deep summary:
    all_analyzed = [analyze_hits(dataset, alb_db, 10000) for dataset in all_hits]
    for fname, analysis in zip(input_filenames, all_analyzed):
        write_analyzed_alb_records(analysis.values(), os.path.join(output_dir, fname + "_analysis.csv"))
        draw_histogram_of_analysis(analysis, os.path.join(output_dir, fname))
        write_analyzed_hits_into_bed_proteome(os.path.join(output_dir, fname + ".proteome.filter_%d.bed" % read_depth_filter), analysis.values())
    
    compare_insertion_and_neighborhood((input_filenames, all_analyzed),
                                       os.path.join(output_dir, "insertion_vs_neighborhood_correlations.txt"))    
    
    # Statistics:
    with open(os.path.join(args.output_dir, "stats.csv"), 'w') as stats_file:
        writer = csv.writer(stats_file, delimiter=',')
        writer.writerow(["File name"] + ALL_STATS)
        for fname, dataset in zip(input_filenames, all_hits):
            stats = get_statistics(dataset, alb_db)
            writer.writerow([fname] + [stats[col] for col in ALL_STATS])
    
    for fname, dataset in zip(input_filenames, all_hits):
        write_hits_into_bed(os.path.join(output_dir, fname + ".filter_%d.bed" % read_depth_filter), dataset)