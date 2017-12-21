# A module based on "ExperimentCoparison.py"
# Hopefully useable in a more universal way.

import os
import argparse
import csv
import glob
import GenomicFeatures
import Shared
from SummaryTable import read_hit_files, analyze_hits, write_data_to_csv

usage = '''AnalyzeAssays.py
    -1st    --first-input-dir   [str]   REQUIRED. Input directory of first *_Hits.txt file.
    -2nd    --second-input-dir  [str]   REQUIRED. Input directory of second *_Hits.txt file.
    -out    --output-dir        [str]   Output directory. Defaults to current directory if left unspecified.
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
        
def compare_two_analyzed_datasets(first_name, first_dataset, second_name, second_dataset, out_dir):
    combined_dataset = []
    
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
    

    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage=usage)
    
    parser.add_argument("-1st", "--first-input-dir", required=True)
    parser.add_argument("-2nd", "--second-input-dir", required=True)
    parser.add_argument("-out", "--output-dir", default='.')
    
    args = parser.parse_args()
    
    first_input_dir = args.first_input_dir
    second_input_dir = args.second_input_dir
    output_dir = args.output_dir
    
    Shared.make_dir(output_dir)
       
    first_name_path = get_filename_path(first_input_dir)
    first_name = first_name_path[0]    
    second_name_path = get_filename_path(second_input_dir)
    second_name = second_name_path[0]

    first_analyzed_tupled = analyze_folder(first_name_path)
    first_analyzed_dataset = first_analyzed_tupled[0]
    second_analyzed_tupled = analyze_folder(second_name_path)
    second_analyzed_dataset = second_analyzed_tupled[0]

    compare_two_analyzed_datasets(first_name, first_analyzed_dataset, second_name, second_analyzed_dataset, output_dir)  
