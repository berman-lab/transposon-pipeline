import os
import csv
import GenomicFeatures
import argparse
import glob
import scipy.stats
import numpy as np
import math
import itertools
import pandas as pd
import Shared
import Organisms
from RangeSet import RangeSet

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle


usage = '''SummaryTable.py  
    -i  --input_dir             [str]   Input directory of '*_Hits.txt'. Defaults to current directory if left unspecified.
    -o  --output-dir            [str]   Output directory for results. Defaults to current directory if left unspecified.
    -f  --read-depth-filter     [str]   Read depth below which insertion events will be ignored. Default is 1
    -h  --help                          Show this help message and exit 
'''

def read_hit_files(files, read_depth_filter=1):
    """Read in the list of hits files.
    
    Parameters
    ----------
    read_depth_filter : int
        The read depth below which insertion events will be ignored.
    """

    return [read_hit_file(f, read_depth_filter) for f in files]

def read_hit_file(filename, read_depth_filter=1):
    """Read in the given hit file.
    
    Parameters
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
            chrom, source, _up_feature_type, up_feature_name, up_gene_name, \
                   up_feature_dist, _down_feature_type, down_feature_name, \
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

def read_pombe_hit_file(filename, read_depth_filter=1):
    # TODO: should be standardized with the Calb hit file.
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

def get_hits_from_wig(wig_file):
    """Read hits from a .wig file.
    
    Used in reading Kornmann's cerevisiae hit data.
    """
    
    cer_db = GenomicFeatures.default_cer_db()
    
    result = []
    with open(wig_file, 'r') as in_file:
        in_file.readline() # Drop the header line
        for line in in_file:
            if line.startswith("variableStep"):
                chrom_name = line[line.index("chrom=") + len("chrom="):].strip()
                continue
             
            if chrom_name == 'chrM':
                continue
            
            hit_pos, reads = line.split()
            hit_pos = int(hit_pos)
            reads = int(reads)
            
            fs = cer_db[chrom_name][hit_pos]
            if len(fs) == 0:
                name = ig_type = "nan"
            else:
                f = fs[0]
                name = f.standard_name
                ig_type = "ORF" # TODO: how do we know it's an ORF? Do we need this?
            
            result.append({"chrom": chrom_name, "hit_pos": hit_pos, "gene_name": name,
                           "ig_type": ig_type,
                           "hit_count": reads})
            
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
    
    alb_db = GenomicFeatures.default_alb_db()
    
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
    
    Parameters
    ----------
    dataset : list of dict
        The hit list, requires that each hit have a "hit_pos", "hit_count" and
        "gene_name" field.
    feature_db : `GenomicFeatures._FeatureDB`
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
        
        # We will count up all of the hits, and then separate them into intra-
        # and inter-feature arrays (for the neighborhood index). The separation
        # is done via using an exon mask:
        exon_mask = np.zeros((chrom_len + 1,), dtype=np.bool)
        for feature in feature_db[chrom]:
            for exon in feature.exons:
                exon_mask[exon[0]:exon[1]+1] = True
                
        domain_mask = np.zeros((chrom_len + 1,), dtype=np.bool)
        for feature in feature_db[chrom]:
            for domain in feature.domains:
                domain_mask[domain[0]:domain[1]+1] = True
        
        # The mask is inverted - the 1 cells are to be kept (they exist and
        # have a high mapq), and the 0 cells are to be ignored.
        # TODO: refactor.
        ignored_mask = np.ones((chrom_len + 1,), dtype=np.bool)
        if isinstance(feature_db, GenomicFeatures.AlbicansFeatureDB):
            ignored_range_set = Organisms.alb.ignored_regions[chrom]
        elif isinstance(feature_db, GenomicFeatures.PombeFeatureDB):
            ignored_range_set = Organisms.pom.ignored_regions[chrom]
        elif isinstance(feature_db, GenomicFeatures.CerevisiaeFeatureDB):
            ignored_range_set = Organisms.cer.ignored_regions[feature_db.get_std_chrom_name(chrom)]
        else:
            ignored_range_set = RangeSet()
        for start, stop in ignored_range_set:
            ignored_mask[start:stop+1] = False
        
        hits_across_chrom = np.zeros((chrom_len + 1,), dtype=np.int)
        reads_across_chrom = np.zeros((chrom_len + 1,), dtype=np.int)
        
        for hit in hits:
            reads_across_chrom[hit["hit_pos"]] += hit["hit_count"]
            # TODO: why do we cap the hits at 2? What happens with pooled
            # hit all_hits?
            hit_pos = hit["hit_pos"]
            if hits_across_chrom[hit_pos] < 2:
                hits_across_chrom[hit_pos] += 1
            
            # TODO: have the gene_name be in standard format - common names can
            # be duplicates, e.g. PDR17
        
        # Hits in ignored locations:
        # There seem to be a few tens of these on each chromosome. A random inspection
        # with IGV suggests that these are bowtie2 alignment errors - it seems to give
        # a high MAPQ to reads that should have a low MAPQ (according to our simulated
        # BAM for finding homologous regions).
        hits_across_chrom = hits_across_chrom * ignored_mask
        reads_across_chrom = reads_across_chrom * ignored_mask
        
        hits_in_features = hits_across_chrom * exon_mask
        hits_outside_features = hits_across_chrom * np.invert(exon_mask)
        reads_outside_features = reads_across_chrom * np.invert(exon_mask)
        
        records = {} # The per-feature records for this chromosome.

        for feature in feature_db[chrom]:
            feature_mask = exon_mask[feature.start:feature.stop+1]
            domain_feature_mask = domain_mask[feature.start:feature.stop+1]
            hits_in_feature = hits_in_features[feature.start:feature.stop+1][feature_mask]
            hits_in_domains_count = hits_in_features[feature.start:feature.stop+1][domain_feature_mask].sum()
            hits_in_feature_count = hits_in_feature.sum()
            reads_in_feature = reads_across_chrom[feature.start:feature.stop+1][feature_mask].sum()
            reads_in_domains = reads_across_chrom[feature.start:feature.stop+1][domain_feature_mask].sum()
            
            # The borders of the neighborhood window:
            window_start = max(1, feature.start - neighborhood_window_size)
            window_end = min(chrom_len, feature.stop + neighborhood_window_size)
            
            hits_outside_feature = hits_outside_features[window_start:window_end+1]
            hits_outside_feature_count = hits_outside_feature.sum()
            
            intergenic_region = feature_db.get_interfeature_range(chrom, (window_start, window_end)) - ignored_range_set
            insertion_index = float(hits_in_feature_count) / feature.coding_length
            
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
                reads_ni = 0
            else:
                neighborhood_insertion_index = float(hits_outside_feature_count) / intergenic_region.coverage
                neighborhood_index = insertion_index / neighborhood_insertion_index
                reads_ni = (
                    (float(reads_in_feature) / feature.coding_length) /
                    (float(reads_outside_features[window_start:window_end+1].sum()) / intergenic_region.coverage)
                )
                
            n_term_insertions = {}
            for n_term_len in (50, 100):
                if feature.strand == 'W':
                    n_term_hits = hits_in_feature[:n_term_len].sum()
                elif feature.strand == 'C':
                    n_term_hits = hits_in_feature[max(0, len(hits_in_feature)-n_term_len):].sum()
                n_term_insertions["n_term_hits_%d" % n_term_len] = n_term_hits
                n_term_insertions["n_term_ni_%d" % n_term_len] = \
                    (float(n_term_hits) / n_term_len) / neighborhood_insertion_index if \
                    hits_outside_feature_count > 0 else 0
            
            # TODO: We assume that the 100 bp upstream is always intergenic, but this isn't always the case.
            # How should we handle this?
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
            upstream_50_hits = upstream_region_50.sum()
            upstream_100_hits = upstream_region_100.sum()
            
            # Compute the longest area without any hits:
            hit_ixes = [0] + list(np.where(hits_in_feature > 0)[0]+1) + [feature.coding_length + 1]
            
            longest_free_intervals = []
            freedom_indices = []
            kornmann_indices = []
            logit_fis = []
            max_skip_tns = 4+1
            for skip_tns in range(max_skip_tns):
                longest_interval = \
                    max(right - left for left, right in
                        zip(hit_ixes, hit_ixes[1+skip_tns:] or [feature.coding_length + 1])) \
                    - 1
                
                longest_free_intervals.append(longest_interval)
                
                freedom_index = float(longest_interval) / feature.coding_length
                freedom_indices.append(freedom_index)
                
                if longest_interval >= 300 and 0.1 < freedom_index < 0.9:
                    kornmann_index = (longest_interval * hits_in_feature_count) / (feature.coding_length ** 1.5)
                else:
                    kornmann_index = 0
                kornmann_indices.append(kornmann_index)
                
                logit_fis.append(freedom_index / (1.0 + math.e ** (-0.01 * (feature.coding_length - 200))))
            
            records[feature.standard_name] = {
                "feature": feature,
                "length": len(feature),
                "hits": hits_in_feature_count,
                "reads": reads_in_feature,
                # Can be tested downstream for zero:
                "neighborhood_hits":  hits_outside_feature_count,
                "nc_window_len":  intergenic_region.coverage,
                "insertion_index": insertion_index,
                "neighborhood_index": neighborhood_index,
                "reads_ni": reads_ni,
                "upstream_hits_50": upstream_50_hits,
                "upstream_hits_100": upstream_100_hits,
                "max_free_region": longest_free_intervals[0],
                "freedom_index": freedom_indices[0],
                "logit_fi": logit_fis[0],
                # The value is singular, to get the coefficient we subtract the pre from the post later on:
                "s_value": log2(reads_in_feature + 1) - total_reads_log, # Add 1 so as to not get a log 0
                # Note: the positions are relative to the gene, and the introns are excised:
                "hit_locations": [ix+1 for (ix, hit) in enumerate(hits_in_feature) if hit > 0],
                "longest_interval": longest_free_intervals[4],
                "kornmann_domain_index": kornmann_indices[4],
                "domain_ratio": float(feature.domains.coverage)  / feature.coding_length,
                "hits_in_domains": hits_in_domains_count,
                "reads_in_domains": reads_in_domains,
                "domain_coverage": feature.domains.coverage,
                "skip_longest_free_intervals": longest_free_intervals,
                "skip_freedom_indices": freedom_indices,
                "skip_kornmann_indices": kornmann_indices,
                "logit_fis": logit_fis[0],
                "bps_between_hits_in_neihgborhood": intergenic_region.coverage / hits_outside_feature_count if hits_outside_feature_count > 0 else 9999,
                "bps_between_hits_in_feature": len(feature) / hits_in_feature_count if hits_in_feature_count > 0 else 9999
            }
            
            records[feature.standard_name].update(n_term_insertions)
            
            for skip_tns in range(max_skip_tns):
                records[feature.standard_name]["longest_free_interval_%d" % skip_tns] = longest_free_intervals[skip_tns]
                records[feature.standard_name]["freedom_index_%d" % skip_tns] = freedom_indices[skip_tns]
                records[feature.standard_name]["kornmann_index_%d" % skip_tns] = kornmann_indices[skip_tns]
                records[feature.standard_name]["logit_fi_%d" % skip_tns] = logit_fis[skip_tns]
            
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
            "field_name": "logit_fi",
            "csv_name": "Logit FI",
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
    
    df = write_data_to_data_frame(data, cols_config)
    df.to_csv(output_filename)

def write_data_to_data_frame(data, cols_config):
    field_col = "field_name"
    format_col = "format"
    csv_name_col = "csv_name"
    sort_by_col = "sort_by"
    
    # Filter out field names that don't exist:
    if data:
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
        # NB: previously, when writing directly to CSV, without the pandas
        # intermediary, we would use `"%s" % s` formatting, but that doesn't
        # work with pandas anymore, so we pass the value as is.
        # We can consider adding a special pandas format field.
        if format_col not in col:
            format_cache[col_name] = lambda s: s
        elif callable(col[format_col]):
            format_cache[col_name] = col[format_col]
        else:
            format_cache[col_name] = lambda s: s
    
    result = pd.DataFrame(
        data=[[format_cache[col[csv_name_col]](record[col[field_col]]) for col in cols_config] for record in ordered_dataset],
        columns=[col[csv_name_col] for col in cols_config]
    )
    
    return result


@Shared.memoized
def get_calb_ess_in_sc():
    cer_essentials = Organisms.cer.literature_essentials
    cer_non_essentials = Organisms.cer.literature_non_essentials
    all_alb_fs = GenomicFeatures.default_alb_db().get_all_features()
    cer_db = GenomicFeatures.default_cer_db()
    
    return (
        set(f.standard_name for f in all_alb_fs if any(cer_db.get_feature_by_name(c).standard_name in cer_essentials for c in f.cerevisiae_orthologs)),
        set(f.standard_name for f in all_alb_fs if any(cer_db.get_feature_by_name(c).standard_name in cer_non_essentials for c in f.cerevisiae_orthologs)),
    )

@Shared.memoized
def get_calb_orths_in_sp():
    return Organisms.get_calb_orths_in_sp()

@Shared.memoized
def get_calb_ess_in_sp():
    # TODO: refactor into the Organisms module.
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
    
    essentials = set()
    non_essentials = set()
    for alb_feature in GenomicFeatures.default_alb_db().get_all_features():
        ortholog_row = joined_table[joined_table["albicans standard name"] == alb_feature.standard_name]
        if ortholog_row.empty:
            continue 
        
        ess_annotation = ortholog_row["essentiality"].iloc[0]
        if ess_annotation == "viable":
            non_essentials.add(alb_feature.standard_name)
        elif ess_annotation == "inviable":
            essentials.add(alb_feature.standard_name)
            
    return (essentials, non_essentials)

def enrich_with_pombe(records):
    """Add pombe orthologs and essentiality annotations to albicans records."""
    
    # TODO: maybe we should add this to analyze_hits? 
    
    essentials, non_essentials = get_calb_ess_in_sp()
    sp_orthologs = get_calb_orths_in_sp()
    
    for record in records:
        alb_name = record["feature"].standard_name
        
        ortholog_name = sp_orthologs.get(alb_name, "")
        ortholog_essentiality = "Yes" if alb_name in essentials else \
                                "No" if alb_name in non_essentials else \
                                ""
         
        record["pombe_ortholog"] = ortholog_name
        record["essential_in_pombe"] = ortholog_essentiality

def enrich_alb_records(records):
    """Enrich albicans records with 3-rd party data.
    
    3-rd party data include orthologs, essentiality from literature, etc."""
    
    # TODO: maybe we should add this to analyze_hits?
    
    cer_db = GenomicFeatures.default_cer_db()
    
    # Pombe:
    enrich_with_pombe(records)

    # Cerevisiae:
    cer_essentials = Organisms.cer.literature_essentials
    cer_non_essentials = Organisms.cer.literature_non_essentials

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
    
    alb_db = GenomicFeatures.default_alb_db()
    alb_pheno_essential_features = set(f.standard_name for f in filter(None, (alb_db.get_feature_by_name(f) for f in alb_pheno_essentials)))
    alb_pheno_non_essential_features = set(f.standard_name for f in  filter(None, (alb_db.get_feature_by_name(f) for f in alb_pheno_non_essentials)))
    
    alb_pheno_toss_up = alb_pheno_essential_features & alb_pheno_non_essential_features
    alb_pheno_essential_consensus = alb_pheno_essential_features - alb_pheno_non_essential_features
    alb_pheno_non_essential_consensus = alb_pheno_non_essential_features - alb_pheno_essential_features
    
    roemer_grace_essentials, roemer_grace_non_essentials, \
        omeara_grace_essentials, omeara_grace_non_essentials = get_grace_essentials()
    
    homann_deletions = get_homann_deletions()
    noble_deletions = get_noble_deletions()
    sanglard_deletions = get_sanglard_deletions()
    
    mitchell_essentials = get_mitchell_essentials()
    
    for record in records:
        feature = record["feature"]
        std_name = feature.standard_name
        if feature.cerevisiae_orthologs:
            # TODO: according to albicans feature file, a gene may have more than one ortholog, but in practice this doesn't happen.
            assert len(feature.cerevisiae_orthologs) == 1
            cer_ortholog = cer_db.get_feature_by_name(iter(feature.cerevisiae_orthologs).next()).standard_name
            
            # TODO: ? should be reserved for conflicting entries, not just those
            # for which data in not available.
            record["essential_in_cerevisiae"] = "Yes" if cer_ortholog in cer_essentials else \
                                                "No" if cer_ortholog in cer_non_essentials else \
                                                "?" if cer_ortholog in Organisms.cer.conflicting_essentials else ""
             
            record["cer_synthetic_lethal"] = cer_synthetic_lethals.get(cer_ortholog, "")
            record["cer_fitness"] = cer_fitness.get(cer_ortholog, "")
        else:
            record["essential_in_cerevisiae"] = ""
            record["cer_synthetic_lethal"] = ""
            record["cer_fitness"] = ""
        record["essential_in_albicans"] = "Yes" if std_name in alb_pheno_essential_consensus else \
            "No" if std_name in alb_pheno_non_essential_consensus else \
            "?" if std_name in alb_pheno_toss_up else ""
        record["essential_in_albicans_grace_roemer"] = "Yes" if std_name in roemer_grace_essentials else \
            "No" if std_name in roemer_grace_non_essentials else ""
        record["essential_in_albicans_grace_omeara"] = "Yes" if feature.standard_name in omeara_grace_essentials else \
            "No" if std_name in omeara_grace_non_essentials else ""
            
        record["homann_deletions"] = "Deleted" if std_name in homann_deletions else ""
        record["noble_deletions"] = "Deleted" if std_name in noble_deletions else ""
        record["sanglard_deletions"] = "Deleted" if std_name in sanglard_deletions else ""
        record["deleted_in_calb"] = ",".join(filter(None, [record["homann_deletions"] and "Homann", record["noble_deletions"] and "Noble", record["sanglard_deletions"] and "Sanglard"]))
        
        record["essential_in_mitchell"] = "Yes" if std_name in mitchell_essentials else ""

@Shared.memoized
def get_homann_deletions():
    # TODO: the deletions should be listed in the Organisms module instead of here.
    alb_db = GenomicFeatures.default_alb_db()
    mutants_filepath = Shared.get_dependency("albicans/MUTANT COLLECTIONS IN PROGRESS July 5 2017.xlsx")
    mutants_table = pd.read_excel(mutants_filepath, skiprows=1, header=None, usecols="A")
    return set(f.standard_name for f in (alb_db.get_feature_by_name(str(n) + "_A") for n in mutants_table[0]) if f)

@Shared.memoized
def get_noble_deletions():
    alb_db = GenomicFeatures.default_alb_db()
    mutants_filepath = Shared.get_dependency("albicans/MUTANT COLLECTIONS IN PROGRESS July 5 2017.xlsx")
    mutants_table = pd.read_excel(mutants_filepath, skiprows=1, header=None, usecols="P")
    return set(f.standard_name for f in (alb_db.get_feature_by_name(n.lower()) for n in mutants_table[0]) if f)    

@Shared.memoized
def get_sanglard_deletions():
    alb_db = GenomicFeatures.default_alb_db()
    mutants_filepath = Shared.get_dependency("albicans/MUTANT COLLECTIONS IN PROGRESS July 5 2017.xlsx")
    mutants_table = pd.read_excel(mutants_filepath, skiprows=1, header=None, usecols="J")
    return set(f.standard_name for f in (alb_db.get_feature_by_name(str(n) + "_A") for n in mutants_table[0]) if f)

@Shared.memoized
def get_grace_essentials():
    """Get names of essential and non-essential feature names in the GRACE
    collection.
    
    Feature names are in the standard format.
    
    Returns
    -------
    tuple of sets of str : (roemer_essentials, roemer_non_essentials,
                            omeara_essentials, omeara_non_essentials)
    """
    
    # TODO: should be listed in the Organisms module instead of here.
    # TODO: figure out if Roemer data actually exists in O'Meara's tables,
    # or was it but an illusion.
    
    grace_table_path = Shared.get_dependency("albicans/ncomms7741-s2.xls")
    grace_table = pd.read_excel(grace_table_path, sheet_name=1, #"Essentiality scores",
                                skiprows=14, header=None,
                                names=["orf19 name", "Common", "Description", "Plate", "Position",
                                       "tet growth", "5-FOA excision", "dox growth", "Essentiality Concordance  (Y/N)",
                                       "Essential/Not essential", "S. cerevisiae homolog", "S. cerevisiae KO phenotype"])
    
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
        
    alb_db = GenomicFeatures.default_alb_db()
    
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
            
    return (roemer_grace_essentials, roemer_grace_non_essentials,
            omeara_grace_essentials, omeara_grace_non_essentials)

@Shared.memoized
def get_mitchell_essentials():
    # TODO: should be in Organisms.
    alb_db = GenomicFeatures.default_alb_db()
    
    with open(Shared.get_dependency("albicans", "List of possibly ess genes from Aaron Mitchell.csv"), 'r') as in_file:
        orfs_19 = [s.strip() for s in in_file.read().split()]
        orfs_22 = filter(None, (alb_db.get_feature_by_name(o) for o in orfs_19))
        
    return set(o.standard_name for o in orfs_22)

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

def perform_pairwise_correlations(analysis_labels, analyzed_records, output_file_prefix):
    """Compute and output Pearson and Spearman correlations between all pairs
    of analyzed hits. Various data points per feature are considered for
    correlation from within the dataset (hits, reads, etc.).
    
    Parameters
    ----------
    analysis_labels : list of strings
        The dataset labels.
    analyzed_records : list of lists
        A list of the datasets, each dataset a list of records.
    output_file_prefix : str
        The prefix of the output files.
    """
    
    zipped_input = sorted(zip(analysis_labels, analyzed_records))
    
    corr_data_sources = ("hits_linear", "hits_log", "reads_linear", "reads_log", "ni_linear", "ni_log")
    corr_types = ("spearman", "pearson")
    
    # Prepare the data to correlated, assume all analyzed records contain the
    # same records.
    corr_data = {corr_data_source: {} for corr_data_source in corr_data_sources}
    for label, records in zipped_input:
        records.sort(key=lambda r: r["feature"].name)
        corr_data["hits_linear"][label] = [r["hits"] for r in records]
        corr_data["reads_linear"][label] = [r["reads"] for r in records]
        corr_data["ni_linear"][label] = [r["neighborhood_index"] for r in records]
        corr_data["hits_log"][label] = [math.log(v) if v > 0 else 0 for v in corr_data["hits_linear"][label]]
        corr_data["reads_log"][label] = [math.log(v) if v > 0 else 0 for v in corr_data["reads_linear"][label]]
        corr_data["ni_log"][label] = [math.log(v) if v > 0 else 0 for v in corr_data["ni_linear"][label]]
        
    # Perform the comparisons:
    corr_results = {corr_data_source: {corr_type: {} for corr_type in corr_types}
                    for corr_data_source in corr_data_sources}
    for (label1, records1), (label2, records2) in itertools.combinations(zipped_input, 2):
        for corr_data_source in corr_data_sources:
            corr_results[corr_data_source]["pearson"][(label1, label2)] = \
                scipy.stats.pearsonr(corr_data[corr_data_source][label1],
                                     corr_data[corr_data_source][label2])
            corr_results[corr_data_source]["spearman"][(label1, label2)] = \
                scipy.stats.spearmanr(corr_data[corr_data_source][label1],
                                      corr_data[corr_data_source][label2])
    
    # Write out everything:
    ordered_labels = [p[0] for p in zipped_input]
    for corr_data_source, corr_type in itertools.product(corr_data_sources, corr_types):
        with open(output_file_prefix + ".%s.%s.csv" % (corr_data_source, corr_type), 'w') as out_file:
            writer = csv.writer(out_file)
            writer.writerow([""] + ordered_labels)
            results = corr_results[corr_data_source][corr_type]
            for ix, row_label in enumerate(ordered_labels[:-1]):
                # NB: not writing out p-values, as most of them are > 0.001,
                # and it's easier to read without them. However, BEWARE! Not
                # all of them are so small!
                writer.writerow([row_label] +
                                ["%.3f" % results[(col_label, row_label)][0] for col_label in ordered_labels[:ix]] +
                                [""] +
                                ["%.3f" % results[(row_label, col_label)][0] for col_label in ordered_labels[ix+1:]])        

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
    alb_db = GenomicFeatures.default_alb_db()
    
    chroms = list(alb_db)
    chroms.sort(key=lambda c: c.name)
    max_chrom_len = max(len(c) for c in chroms)
    
    fig = plt.figure(figsize=(30, 4*len(chroms)))
    plt.gcf().suptitle(title, fontsize=22, fontweight='bold', y=0.92, x=0.2)
    
    ade2 = alb_db.get_feature_by_name("ADE2")
    
    # How many chromosomes should we display?
    row_num = len(chroms)
    for chrom in chroms:
        d = data[chrom.name]
        if max(d) > y_max:
            row_num += 1
    
    # Start drawing each chromosome:
    base_row = 0
    for chrom in sorted(chroms, key=lambda c: int(c.name[7]) if c.name[7] != "R" else 0):
        chrom_display_name = chrom.name[4].capitalize() + chrom.name[5:8]
        chrom_data = data[chrom.name]
        
        if max(chrom_data) <= y_max:
            # The data in the chromosome does not exceed the y maximum, draw
            # it as usual.
            ax = plt.subplot2grid((row_num, 100),
                                  (base_row, 0),
                                  colspan=len(chrom)*100/max_chrom_len)
            ax.set_ylim([0, y_max])
            base_row += 1
        elif max(chrom_data) <= y_max * 2:
            # The data in the chromosome is less than twice as large as the
            # maximum, allocate it twice the space and draw as usual. 
            ax = plt.subplot2grid((row_num, 100),
                                  (base_row, 0),
                                  rowspan=2,
                                  colspan=len(chrom)*100/max_chrom_len)
            ax.set_ylim([0, y_max*2])
            base_row += 2
        else:
            # The data in the chromosome is too large (twice than the maximum y
            # value), and so will undergo a jump-cut along the y axis.
            ax = plt.subplot2grid((row_num, 100),
                                  (base_row+1, 0),
                                  colspan=len(chrom)*100/max_chrom_len)
            ax.set_ylim([0, y_max])
            ax2 = plt.subplot2grid((row_num, 100),
                                   (base_row, 0),
                                   colspan=len(chrom)*100/max_chrom_len,
                                   sharex=ax)
            max_top_y = max(chrom_data) / 100 * 100 + 100
            top_y = [max_top_y - y_max, max_top_y]
            base_row += 2
            
            ax2.set_ylim(top_y)
            ax2.bar(range(len(chrom_data)), chrom_data)
            
            ax2.spines['bottom'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax2.xaxis.tick_top()
            ax2.tick_params(labeltop='off')  # don't put tick labels at the top
            ax.xaxis.tick_bottom()
            
            d = .015  # how big to make the diagonal lines in axes coordinates
            # arguments to pass to plot, just so we don't keep repeating them
            kwargs = dict(transform=ax2.transAxes, color='k', clip_on=False)
            ax2.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
            ax2.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
            
            kwargs.update(transform=ax.transAxes)  # switch to the bottom axes
            ax.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
            ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
        
        # Adding the last position as 0/1 so the plot won't be shorter:    
        ax.bar(range(len(chrom_data)), chrom_data)
        ax.set_ylabel(chrom_display_name, fontsize = 20, fontweight = 'bold')
        ax.set_xlim([0, len(chrom_data)])
        ax.tick_params(axis='both', which='major', labelsize=14)
        
        patches = {
            "centromere": "red",
            "repeat_region": "green",
            "long_terminal_repeat": "cyan",
            "tRNA|Verified": "orange",
            "tRNA|Uncharacterized": "orange",
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
                feature_len = feature.coding_length / 3
                aa_hits = [feature_len - h + 1 for h in aa_hits]
            for aa_hit in aa_hits:
                bed_file.write("%s\t%d\t%d\n" % (record["feature"].standard_name, aa_hit, aa_hit+1))

def make_hit_map(hits, bin_size):
    alb_db = GenomicFeatures.default_alb_db()
    
    result = {chrom.name: [0]*(len(chrom)/bin_size + 1) for chrom in alb_db}
    
    for hit in hits:
        result[hit["chrom"]][hit["hit_pos"]/bin_size] += 1
    
    return result

def make_read_map(hits, bin_size):
    alb_db = GenomicFeatures.default_alb_db()
    
    result = {chrom.name: [0]*(len(chrom)/bin_size + 1) for chrom in alb_db}
    
    for hit in hits:
        result[hit["chrom"]][hit["hit_pos"]/bin_size] += hit["hit_count"]
    
    return result

def transform_chrom_map_with_log(read_map, base):
    for data in read_map.values():
        for ix in range(len(data)):
            data[ix] = math.log(data[ix], base) if data[ix] > 0 else 0
    
    return read_map

if __name__ == "__main__":    
    parser = argparse.ArgumentParser(usage=usage)
    
    parser.add_argument("-i", "--input-dir", default='.')
    parser.add_argument("-o", "--output-dir", default='.')
    parser.add_argument("-f", "--read-depth-filter", type=int, default=1)
    args = parser.parse_args()
    
    input_dir = args.input_dir
    output_dir = args.output_dir
    read_depth_filter = args.read_depth_filter
    
    Shared.make_dir(output_dir)
    
    alb_db = GenomicFeatures.default_alb_db()
    
    input_file_paths = glob.glob(os.path.join(input_dir, "*_Hits.txt"))
    input_filenames = [os.path.split(file_path)[-1][:-9] for file_path in input_file_paths]
    
    all_hits = read_hit_files(input_file_paths, read_depth_filter)
    
    # Write per-bin hits and reads, for 3D analysis:
    bin_size = 10000
    # [hits, reads, hit rank, read rank]
    hit_read_bins = {fname: {chrom.name: [[0, 0, 0, 0] for _ in range(len(chrom) / bin_size + 1)] for chrom in alb_db} for fname in input_filenames}
    for fname, hits in zip(input_filenames, all_hits):
        fname_hit_read_bins = hit_read_bins[fname]
        for hit in hits:
            bin = hit["hit_pos"] / bin_size
            bin_data = fname_hit_read_bins[hit["chrom"]][bin]
            bin_data[0] += 1
            bin_data[1] += hit["hit_count"]
    
    from itertools import chain
    for fname_hit_read_bins in hit_read_bins.values():
        for bin_rank, bin_data in enumerate(sorted(chain(*fname_hit_read_bins.values()), key=lambda b: b[0], reverse=True)):
            bin_data[2] = bin_rank+1
            
    with open(os.path.join(output_dir, "binned_hits.RDF_%d.csv" % read_depth_filter), "w") as out_file:
        writer = csv.writer(out_file)
        writer.writerow(["Bin index"] + list(chain(*zip(input_filenames, ["%s rank" % fname for fname in input_filenames]))))
        import re
        for chrom in sorted(alb_db, key=lambda c: c.name):
            chrom_ix = re.match(".*chr(.)", chrom.name).group(1)
            for bin_ix in range(len(chrom) / bin_size + 1):
                writer.writerow(["%s-%d" % (chrom_ix, bin_ix)] + list(chain(*[hit_read_bins[fname][chrom.name][bin_ix][0::2] for fname in input_filenames])))
    
    all_analyzed = [analyze_hits(dataset, alb_db, 10000) for dataset in all_hits]
    
    # Write nominal hit-per-ORF table (one table for all datasets):
    with open(os.path.join(output_dir, "hit_summary.RDF_%d.csv" % read_depth_filter), 'w') as out_file:
        writer = csv.writer(out_file)
        writer.writerow(["Standard name", "Common name"] + [f[-2:] for f in input_filenames])
        for record_row in zip(*[ds.values() for ds in all_analyzed]):
            assert len(set(r["feature"].standard_name for r in record_row)) == 1
            feature = record_row[0]["feature"]
            writer.writerow([feature.standard_name, feature.common_name] +
                            [r["hits"] for r in record_row])
    
    # Write a table of reads-per-hit for each dataset:
    for fname, hits in zip(input_filenames, all_hits):
        # Transform:
        sum_hits = {chrom.name: {} for chrom in alb_db}
        for h in hits:
            sum_hits[h["chrom"]][h["hit_pos"]] = sum_hits[h["chrom"]].get(h["hit_pos"], 0) + h["hit_count"] 
            
        # Write out:
        with open(os.path.join(output_dir, "all_hits.%s.csv" % fname), 'w') as out_file:
            writer = csv.writer(out_file)
            writer.writerow(["Chromosome", "Position", "Reads"])
            for chrom, data in sorted(sum_hits.items()):
                for pos, reads in sorted(data.items()):
                    writer.writerow([chrom, pos, reads])
    
    # Plot the read distributions per hit:
    for fname, hits in zip(input_filenames, all_hits):
        plt.hist([math.log(h["hit_count"]) for h in hits], bins=100, range=(0,10))
        plt.title("Read distribution per hit in %s" % fname)
        plt.xlabel("log10 of # of reads per hit")
        plt.ylabel("# of hits")
        plt.savefig(os.path.join(output_dir, "reads_distribution.hits.log10.%s.rdf_%d.png" % (fname, read_depth_filter)))
        plt.close()
    
    # Plot the read distributions per gene:
    for fname, analysis in zip(input_filenames, all_analyzed):
        plt.hist([math.log(r["reads"] + 1, 10) if r["reads"] > 0 else -0.1 for r in analysis.values()],
                 bins=[-0.1] + [i/10.0 for i in range(61)])
        plt.title("Read distribution per feature in %s" % fname)
        plt.xlabel("log10 of # of reads per feature")
        plt.ylabel("# of features")
        plt.savefig(os.path.join(output_dir, "reads_distribution.log10.%s.rdf_%d.png" % (fname, read_depth_filter)))
        plt.close()
    
    # Plot the hit and read chromosomal maps:
    bin_size = 10000
    for fname, hits in zip(input_filenames, all_hits):
        hit_map = make_hit_map(hits, bin_size)
        read_map = make_read_map(hits, bin_size)
        # NB: using a different log base doesn't affect the final figure much,
        # as the relationships between the values stay the same no matter what
        # base is used. Only the axis labels change.
        log_read_map = transform_chrom_map_with_log(make_read_map(hits, bin_size), 10)
           
        all_hits_sorted = sorted(sum(hit_map.values(), []))
        max_hits = int( round(max(all_hits_sorted), -2) + 100 )
        bottom_cut = int( round(all_hits_sorted[int(len(all_hits_sorted) * 0.5)], -2) + 100 )
        bottom_cut = min(map(max, hit_map.values())) / 100 * 100 + 100
        draw_chrom_map(hit_map, bin_size, "Hits/10000 bps", bottom_cut, os.path.join(output_dir, "hit_map.%s.png" % fname))
          
        draw_chrom_map(log_read_map, bin_size, "Log10 reads/10000 bps", 7, os.path.join(output_dir, "log10_read_map.%s.png" % fname))
        draw_chrom_map(read_map, bin_size, "Reads/10000 bps", 100000, os.path.join(output_dir, "read_map.%s.png" % fname))
    
    # More correlations than you can shake a stick at:
    perform_pairwise_correlations(input_filenames,
                                  [a.values() for a in all_analyzed],
                                  os.path.join(output_dir, "correlations"))
    
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