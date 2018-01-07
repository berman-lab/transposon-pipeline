# A program to compare and analyze two assays
# Should consider moving Pearson and Spearman correlation module from SummaryTable to here
# (This would run all single-file analyses via SummaryTable, and all two-file analyses here)

import os
import argparse
import csv
import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import GenomicFeatures
import Shared
from SummaryTable import read_hit_files, analyze_hits, write_data_to_csv
from matplotlib_venn import venn2

usage = '''AnalyzeAssays.py
    -1st    --first-input-dir   [str]   REQUIRED. Input directory of first *_Hits.txt file.
    -2nd    --second-input-dir  [str]   REQUIRED. Input directory of second *_Hits.txt file.
    -out    --output-dir        [str]   Output directory. Defaults to current directory if left unspecified.
    -d      --histo-dist                Draw histogram of S score distribution.
    -v      --venn-draw                 Draw Venn diagram of C albicans genes with hits.
    -h      --help                      Show this help message and exit 
'''

RDF = 1     # read depth filter

def main():
    parser = argparse.ArgumentParser(usage=usage)
    
    parser.add_argument("-1st", "--first-input-dir", required=True)
    parser.add_argument("-2nd", "--second-input-dir", required=True)
    parser.add_argument("-out", "--output-dir", default='.')
    parser.add_argument("-d", "--histo-dist", default=False, action='store_true')
    parser.add_argument("-v", "--venn-draw", default=False, action='store_true')
    
    args = parser.parse_args()
    
    first_input_dir = args.first_input_dir
    second_input_dir = args.second_input_dir
    output_dir = args.output_dir
    histo_dist = args.histo_dist
    venn_draw = args.venn_draw
    correlations = args.correlations
    
    Shared.make_dir(output_dir)
       
    first_name_path = get_filename_path(first_input_dir)
    first_name = first_name_path[0]    
    second_name_path = get_filename_path(second_input_dir)
    second_name = second_name_path[0]

    first_analyzed_tupled = analyze_folder(first_name_path)
    first_analyzed_dataset = first_analyzed_tupled[0]
    second_analyzed_tupled = analyze_folder(second_name_path)
    second_analyzed_dataset = second_analyzed_tupled[0]

    data = compare_two_analyzed_datasets(first_name, first_analyzed_dataset, second_name, second_analyzed_dataset, output_dir, venn_draw, correlations)

    if venn_draw:
        draw_venn_diag(data, first_name, second_name, output_dir)

    if histo_dist:
        plot_s_scores(first_name, first_analyzed_dataset, second_name, second_analyzed_dataset, output_dir)

def get_filename_path(folder):
    """Gets path and filename for hits file

    Parameters
    ----------
    folder  :   string
        An directory path as input when running the program

    Returns
    -------
    list
        filename for hits file as first element, 
        filepath for hits file as second element
    """
    file_paths = glob.glob(os.path.join(folder, "*_Hits.txt"))
    filenames = [os.path.split(file_path)[-1][:-22] for file_path in file_paths]

    file_path = file_paths[0]
    filename = filenames[0]

    return filename, file_path

def analyze_folder(file_name_path):
    """Runs analyze_hits from SummaryTable using output from get_filename_path

    Parameters
    ----------
    file_name_path  :  list 
        The hits file filename and filepath as elements 0 and 1 respectively 

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
    file_path = [file_name_path[1]]
    filenames = [file_name_path[0]]

    alb_db = GenomicFeatures.default_alb_db()
    all_hits = read_hit_files(file_path, RDF)
    return [analyze_hits(hits, alb_db, 10000).values() 
            for hits in all_hits]
        
def compare_two_analyzed_datasets(first_name, first_dataset, second_name, second_dataset, out_dir, venn_draw): #!!!, correlations):
    """Creates and writes combined record for two input datasets, and optionally gets hit data for Venn diagram

    Parameters
    ----------
    *_name  : string
        Filename for first or second input hits file directory
    *_dataset   :   dict of dicts
        Map of standard feature names to analyzed results for first or second input
    out_dir :   string
        Name for output directory
    venn_draw   :   boolean
        Setting for whether or not to draw Venn diagram.
    correlate   :   boolean
        Setting for whether or not to find correlation data

    Writes
    ------
    s_score_analysi.[first_name]_vs_[second_name].csv   :   csv file
        Combined record with statistical analysis for two inputs in form of .csv file.

    Returns (optional)
    -------
    list of integers
        first_only  :   first element of list
            Number of genes with hits in the first dataset that were not hit in the second
        second_only :   second element
            Number of genes with hits in the second dataset that were not hit in the first
        intersection    :   third element
            Number of genes with hits in both datassets.
    """
    first_only = 0
    second_only = 0
    intersection = 0

    combined_dataset = [first_only, second_only, intersection]
    
    for r1, r2 in zip(first_dataset, second_dataset):
        assert r1['feature'] == r2['feature']   # alternately could use r1.sort(key=lambda r: r["feature"].name), r2.sort(etc.)
        
        combined_record = dict(r1)
        combined_dataset.append(combined_record)
        
        combined_record["s_value_first"] = r1["s_value"]
        combined_record["reads_first"] = r1["reads"]
        combined_record["hits_first"] = r1["hits"]
                
        combined_record["s_value_second"] = r2["s_value"]
        combined_record["reads_second"] = r2["reads"]
        combined_record["hits_second"] = r2["hits"]
        
        combined_record["s_score"] = r2["s_value"] - r1["s_value"]

        if venn_draw:
            if r1['hits'] > 0 and r2['hits'] == 0:
                first_only += 1 
            elif r1['hits'] == 0 and r2['hits'] > 0:
                second_only += 1 
            elif r1['hits'] > 0 and r2['hits'] > 0:
                intersection +=1 
          
    cols_config = [
        {
            "field_name": "feature",
            "csv_name": "Standard name",
            "format": lambda f:f.standard_name
        },
        
        {
            "field_name": "feature",
            "csv_name": "Common name",
            "format": lambda f:f.common_name
        },
        
        {
            "field_name": "hits_first",
            "csv_name": "Hits in %s" % first_name,
            "format": "%d"
        },
        
        {
            "field_name": "hits_second",
            "csv_name": "Hits in %s" % second_name,
            "format": "%d"
        },
        
        {
            "field_name": "feature",
            "csv_name": "Length",
            "format": lambda f: len(f)
        },
        
        {
            "field_name": "reads_first",
            "csv_name": "Reads in %s" % first_name,
            "format": "%d"
        },
        
        {
            "field_name": "reads_second",
            "csv_name": "Reads in %s" % second_name,
            "format": "%d"
        },
        
        {
            "field_name": "s_value_first",
            "csv_name": "S value in %s" % first_name,
            "format": "%.2f"
        },
        
        {
            "field_name": "s_value_second",
            "csv_name": "S value in %s" % second_name,
            "format": "%.2f"
        },
        
        {
            "field_name": "s_score",
            "csv_name": "S score",
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
    
    out_name = "s_score_analysis.%s_vs_%s.csv" % (first_name, second_name)
 
    write_data_to_csv(combined_dataset, cols_config, os.path.join(out_dir, out_name))
    if venn_draw: return first_only, second_only, intersection

def draw_venn_diag(venn_dataset, first_name, second_name, out_dir):
    """Draws Venn diagram comparing genes hit in each dataset.

    Parameters
    ----------
    venn_dataset    :   list of integers
        First element   :   Number of genes with hits in first dataset only
        Second element  :   Number of genes with hits in second dataset only
        Third element   :   Number of genes with hits in both datasets
    *_name  :   string
        Filename for first or second input hits file directory
    out_dir :   string
        Name for output directory

    Writes
    -------
    venn_diag.*.png :   png file
        Venn diagram image
    """
    venn2(subsets=(venn_dataset[0], venn_dataset[1], venn_dataset[2]), 
            set_labels = (first_name, second_name))
    
    plt.title('$C.$ $albicans$ genes with hits')
    plt.savefig(os.path.join(out_dir, 'venn_diag.%s_and_%s.png' % (first_name, second_name)))
    plt.close()
    
def plot_s_scores(first_name, first_dataset, second_name, second_dataset, out_dir):
    """Draws histogram of S score distribution

    Parameters
    ----------
    *_name  : string
        Filename for first or second input hits file directory
    *_dataset   :   dict of dicts
        Map of standard feature names to analyzed results for first or second input
    out_dir :   string
        Name for output directory

    Writes
    -------
    s_score_dist.*.png  :   png file
        Histogram of S score (difference between S values) distribution
    """
    histogram_values = [r2['s_value'] - r1['s_value'] 
                    for r1, r2 in zip(first_dataset, second_dataset)
                        if r1['reads'] > 0 and r2['reads'] > 0]    

    plt.hist(histogram_values, bins=5*16, range=(-8,8))
    plt.title('S score distributions in\n%s vs %s' % (first_name, second_name))
    plt.xlabel('S scores ($\Delta$ of log$_2$(reads) in a feature)')
    plt.ylabel('Number of features')

    plt.savefig(os.path.join(out_dir, 's_score_dist.%s_vs_%s.png' % (first_name, second_name)))
    plt.close()
    
if __name__ == '__main__':
    main()
