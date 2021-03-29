"""
In order for this script to run it is necessary for VEP command line tool
(https://www.ensembl.org/info/docs/tools/vep/script/index.html) to be installed along with full cache file: 
homo_sapiens_vep_102_GRCh38.tar.gz (ftp://ftp.ensembl.org/pub/current_variation/indexed_vep_cache/)

For the pathogenic file an alignment was needed in order to run the VEP command line tool.

ZoomVar script requires a VEP file so both files must be given to rest_api.py as input.

VEP commands
./vep --cache --format id --symbol --canonical --protein -i ~/Desktop/SNPs-Prediction/04_Features_computation/proteins_alignment/refseqs.txt -o non_pathogenic_vep.txt --force
./vep --cache --format id --symbol --canonical --protein -i ~/Desktop/SNPs-Prediction/04_Features_computation/proteins_alignment/pathogenic_refseqs.txt -o pathogenic_vep.txt --force

ZoomVar commands
python3 rest_script.py ~/Desktop/ensembl-vep/non_pathogenic_vep.txt non_pathogenic_structure.csv -v -s -n 20 -i 20
python3 rest_script.py ~/Desktop/ensembl-vep/pathogenic_vep.txt pathogenic_structure.csv -v -s -n 20 -i 20
"""


'''
Necessary files: refseqs.txt and pathogenic_refseqs.txt
'''

import pandas as pd
import os
import re

script_path = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))
parent_dir_path = os.path.abspath(os.path.join(script_path, os.pardir))
features_dir_path = os.path.abspath(os.path.join(parent_dir_path, os.pardir))


def mergeZoomVarFeatures(dataset_file, features_struct_file, output, path_type):
    """Merges a file which contains the computed features, i.e. 11 columns to the previous complete dataset"""

    print(path_type)
    print("Aligning Q(sasa) and struct_annotations features with current dataset!")

    dataset = pd.read_csv(dataset_file, sep='\t',low_memory=False)
    labels = dataset["labels"]

    if path_type == "non pathogenic":
        protein_id = "id" #indicating the id to be compared, RefSeq of Uniprot

    elif path_type == "pathogenic":
        protein_id = "uniprot"


    features_struct = pd.read_csv(features_struct_file, sep=",")[[protein_id, "uniprot_pos",
                                                            "Q(sasa)", "struct_interactions"]]
    features_names = ["Q(sasa)", "struct_interactions"]

    for i in range(2):
        feature_dict = {}  # values will be replaced with the computed values, else they remain None if not found
        for label in labels: feature_dict[label] = "None"  # dictionary initialization
        feature_name = features_names[i]

        for key in feature_dict.keys():
            for feat in features_struct[[protein_id, "uniprot_pos", feature_name]].values:

                protein = feat[0]
                position = feat[1]
                # if protein id and position is found
                if key[:6].split("|")[0] == protein and \
                        int(re.sub("[^0-9]", "", key.split("|")[-1])) == int(position):
                    # when find the same protein id and position we assign the computed values
                    if feature_name == "struct_interactions":
                        if feature_dict[key] == "None": # not assigned yet
                            if feat[2] == "None":
                                feature_dict[key] = "0"

                            else: # if feat[2] contains more than "None"
                                feature_dict[key] = str(len(feat[2].split(",")))

                        else:
                            # if not "None" in feat[2] and the new record is bigger than the previous one
                            if feat[2] != "None" and int(feature_dict[key]) < int(len(feat[2].split(","))):
                                feature_dict[key] = str(len(feat[2].split(",")))

                    else:
                        if feat[2] == "None":
                            feature_dict[key] = "0"
                        else: # has number
                            feature_dict[key] = str(feat[2])

        print(feature_name + " done!")
        dataset[feature_name] = feature_dict.values()
        print(dataset[feature_name].value_counts())
        out = dataset
        out.to_csv(output, sep="\t", index=False)

    print("Q(sasa) and struct_interactions merged with "+path_type+ " dataset!")


if __name__=="__main__":

    print("All files are already created and stored in ZoomVar_query_script/ folder \n ")
    while True:
        ans = input("Do you want to merge files with dataset? y/n")
        if ans == "y" or ans == "Y":

            merge_ans = input("merge pathogenic?")
            if merge_ans == "y" or merge_ans == "Y":

                mergeZoomVarFeatures(dataset_file= features_dir_path + "/pathogenic_dataset.tsv",
                                    features_struct_file = script_path + "/pathogenic_structure.csv",
                                    output= features_dir_path + "/pathogenic_dataset.tsv",
                                    path_type="pathogenic")

            merge_ans = input("merge non pathogenic?")
            if merge_ans == "y" or merge_ans == "Y":

                mergeZoomVarFeatures(dataset_file= features_dir_path + "/non_pathogenic_dataset.tsv",
                                     features_struct_file =script_path + "/non_pathogenic_structure.csv",
                                     output=features_dir_path + "/non_pathogenic_dataset.tsv",
                                     path_type="non pathogenic")

        elif ans == "n" or "N":
            print("Exiting...")
            break
        else:
            print("Wrong inout given, try again by typing y or n.")
