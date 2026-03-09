#!/usr/bin/env python3

import argparse
import pandas as pd
import os
import subprocess

def main():

    url_bige = 'https://raw.githubusercontent.com/ArnaudBelcour/bigecyhmm/refs/heads/main/bigecyhmm/hmm_databases/hmm_table_template.tsv'
    path_t = get_table(output_dir=".",url=url_bige)

    list_kos = parse_table(path_t,category="Carbon fixation")

def get_table(output_dir:str,url:str):

    command_args = [
        "wget", 
        url, 
        "-P", 
        output_dir,
        "-q"
    ]
    filename = url.split('/')[-1]
    file_path = os.path.join(output_dir,filename)

    if os.path.exists(file_path):
        return file_path
    
    else:
    
        try:
            subprocess.run(command_args, check=True)

            table_path = os.path.join(output_dir, filename)

            print(f"Successfully downloaded table to: {table_path}")
            print(table_path.head())
            return table_path

        except Exception as e:
            print(f"Error during download. Check if 'wget' is installed and the URL is correct.")
            return None



def parse_table(table_file:str,category:str)->list:
    
    all_kos = []

    table_ = pd.read_csv(table_file,sep="\t")
    table_cat_ = table_[table_ ["Category"] == category]

    for index,row in table_cat_.iterrows():
        ko_value = row["Corresponding KO"]
        if isinstance(ko_value, str):
            kos = [k.strip() for k in ko_value.split(',')]
        else:
            kos = [str(ko_value)]
            if kos[0] == 'nan':
                continue
        
        all_kos.extend(kos)
    unique_kos = set(all_kos)
    print(len(unique_kos))
    print(unique_kos)

    return unique_kos




def parse_args():
    parser = argparse.ArgumentParser("TEXT HERE")
    
    parser.add_argument("-i", "--metadata_path",
        help="Excel file containing the metadata information.",
        type=str
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
 
