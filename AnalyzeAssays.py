# A module based on "ExperimentCoparison.py"
# Hopefully useable in a more universal way.
# Need help with ven diagram stuff

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
    -d      --histo-dist        [bool]  Draw histogram of S score distribution. Default False.
    -v      --venn-draw         [bool]  Draw venn diagram of C albicans genes with hits. Default False.
    -h      --help                      Show this help message and exit 
'''

RDF = 1 # read depth filter

def get_filename_path(folder):
    file_paths = glob.glob(os.path.join(folder, "*_Hits.txt"))
    filenames = [os.path.split(file_path)[-1][:-22] for file_path in file_paths]

    file_path = file_paths[0]
    filename = filenames[0]

    return filename, file_path

def analyze_folder(file_name_path):
    file_path = [file_name_path[1]]
    filenames = [file_name_path[0]]

    alb_db = GenomicFeatures.default_alb_db()
    all_hits = read_hit_files(file_path, RDF)
    return [analyze_hits(hits, alb_db, 10000).values() 
            for hits in all_hits]
        
def compare_two_analyzed_datasets(first_name, first_dataset, second_name, second_dataset, out_dir, venn_draw):
    combined_dataset = []

    first_only = 0
    second_only = 0
    intersection = 0
    
    for r1, r2 in zip(first_dataset, second_dataset):
        assert r1['feature'] == r2['feature']
        
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
    venn2(subsets=(venn_dataset[0], venn_dataset[1], venn_dataset[2]), 
            set_labels = (first_name, second_name))
    
    plt.title('$C.$ $albicans$ genes with hits')
    plt.savefig(os.path.join(out_dir, 'venn_diag.%s_and_%s.png' % (first_name, second_name)))
    plt.close()
    
def plot_s_scores(first_name, first_dataset, second_name, second_dataset, out_dir):
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
    
    Shared.make_dir(output_dir)
       
    first_name_path = get_filename_path(first_input_dir)
    first_name = first_name_path[0]    
    second_name_path = get_filename_path(second_input_dir)
    second_name = second_name_path[0]

    first_analyzed_tupled = analyze_folder(first_name_path)
    first_analyzed_dataset = first_analyzed_tupled[0]
    second_analyzed_tupled = analyze_folder(second_name_path)
    second_analyzed_dataset = second_analyzed_tupled[0]

    if venn_draw:
        venn_data = compare_two_analyzed_datasets(first_name, first_analyzed_dataset, second_name, second_analyzed_dataset, 
            output_dir, venn_draw)    
        draw_venn_diag(venn_data, first_name, second_name, output_dir)

    else:
        compare_two_analyzed_datasets(first_name, first_analyzed_dataset, second_name, second_analyzed_dataset, 
            output_dir, venn_draw)            

    if histo_dist:
        plot_s_scores(first_name, first_analyzed_dataset, second_name, second_analyzed_dataset, output_dir)
