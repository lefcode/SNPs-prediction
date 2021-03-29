'''
Necessary files: PTMD.txt and whole gzipped_files/ folder
http://dbptm.mbc.nctu.edu.tw/download.php
http://ptmd.biocuckoo.org/download.php
'''

import os
import pandas as pd
from os import listdir
from os.path import isfile, join
import re

script_path = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))
parent_dir_path = os.path.abspath(os.path.join(script_path, os.pardir))
features_dir_path = os.path.abspath(os.path.join(parent_dir_path, os.pardir))


def filterPTMD(ptmd_file):
    """Filters all non-relevant columns and species rather than homo sapiens"""

    ptmd = pd.read_csv(ptmd_file, sep="\t",
                       names=["uni", "gene", "disease", "ptm", "wt", "mut", "pos", "publ", "species"])

    ptmd = ptmd.drop(["gene", "disease", "publ"], axis=1)
    ptmd.drop(ptmd[ptmd['species'] != "Homo sapiens (Human)"].index, inplace=True)
    ptms = ptmd['ptm'].value_counts()

    print(ptms)


def createFeatureFile(list_of_file_names):
    """ Merges all data files into one by first filtering out all non Human species
    It creates ptms_feature_file.tsv with columns:Uniprot  position  PTM """

    with open(script_path + "/ptms_feature_file.tsv", "w+") as features_file:
        features_file.write("Uniprot\tposition\tPTM\n")
        for file in list_of_file_names:
            ptm_data = pd.read_csv(script_path + "/data/" + file, sep="\t",
                                   names=["Species", "Uniprot", "Position", "PTM",
                                          "Reference", "Sequence"]) \
                .drop(["Reference", "Sequence"], axis=1)

            # filter out all the non human records
            for record in ptm_data.values:
                if 'HUMAN' in record[0]:
                    features_file.write(record[1] + "\t" + str(record[2]) + "\t" + record[3] + "\n")

    with open(script_path+"/filtered_PTMD.tsv","w+") as filtered_PTMD:
        # add data from PTMD.txt file
        ptmd = pd.read_csv(script_path + "/PTMD.txt", sep="\t")
        for record in ptmd.values:
            if record[8] == "Homo sapiens (Human)":
                filtered_PTMD.write(record[0] + "\t" + str(record[6]) + "\t" + record[3] + "\n")

    print("Pathogenic feature file created!")
    print("ptms_feature_file.tsv file created!")


def createAlignmentFile(feature_file, filtered_ptmd_file, alignment_file):
    """ With respect to the alignment file that is already created, we create an alignment file
    which contains the same data but with refseq ids """

    feature = pd.read_csv(feature_file, sep="\t", dtype="str")
    filtered_feature = pd.read_csv(filtered_ptmd_file, sep="\t")
    alignment = pd.read_csv(alignment_file, sep=",", names=["refseq", "uniprot"])
    count = 0
    print("creating non pathogenic feature file")

    with open(script_path + "/refseqs_ptms_feature_file.tsv", "w+") as align:
        align.write("refseq\tposition\tPTM\n")
        for feat in feature.values:
            for al in alignment.values:
                if feat[0] == al[1]:  # same uniprot id
                    align.write(al[0] + "\t" + str(feat[1]) + "\t" + feat[2] + "\n")
                    count += 1
                    break


   #with open(script_path+"/filtered_refseqs_PTMD.tsv","w+") as filtered_refseqs_PTMD:
        # filtered_refseqs_PTMD.write("refseq\tposition\tPTM\n")

        for feat in filtered_feature.values:
            for al in alignment.values:
                if feat[0] == al[1]:  # same uniprot id
                    align.write(al[0] + "\t" + str(feat[1]) + "\t" + feat[2] + "\n")
                    count += 1
                    break

    print("Non pathogenic feature file created!")
    return script_path + "/refseqs_ptms_feature_file.tsv"


def createFeatureWithDataset(list_of_freq_ptms, dataset_file, feature_file, output, path_type):
    """ Aligns proteins with PTMs and creates a feature file with 1s and 0s assigned to each protein"""

    print(path_type)
    print("Creating feature data file! It should latter be merged with the dataset without any "
          "processing and filtering as it is completely aligned with it! ")
    dataset = pd.read_csv(dataset_file, sep='\t', low_memory=False)
    feature = pd.read_csv(feature_file, sep="\t", dtype="str")

    labels = dataset["labels"]

    feature_dict = {}  # values will be replaced with the computed values, else they remain None if not found
    for label in labels: feature_dict[label] = "None"  # dictionary initialization

    with open(output, "w+") as ptms_file:
        for i in range(len(list_of_freq_ptms)):
            ptms_file.write(list_of_freq_ptms[i]+ "\t")
        ptms_file.write("Non-Freq-PTM\n")

        for key in feature_dict.keys():
            ptm_flag = 0
            for feat in feature.values:
                protein = feat[0]
                position = feat[1]
                if ptm_flag == 0 and type(position) != float:  # only if not found yet and position!=nan

                    if key.split("|")[0] == protein:
                        positions = position.split('/') # some positions are lists with delimiter "\"
                        for pos in positions:
                            if int(re.sub("[^0-9]", "", key.split("|")[-1])) == int(pos):
                              #  int(re.sub("[^0-9]", "", key.split("|")[1])) == int(position):

                                ptm_flag = 1
                                # if same protein and position found then put 1 to the equivalent PTM and 0 to all
                                # the others
                                freq_ptm_flag = 0

                                for i in range(len(list_of_freq_ptms)):
                                    feature_name = list_of_freq_ptms[i]
                                    if feat[2] == feature_name:
                                        ptms_file.write("1\t")
                                        freq_ptm_flag = 1
                                    else:
                                        ptms_file.write("0\t")


                                # if already found in at least one frequent ptm then assign 0 else 1
                                if freq_ptm_flag == 0:

                                    ptms_file.write("1\n")
                                else:
                                    ptms_file.write("0\n")

            if ptm_flag == 0:  # then it wasn't found fill in with "None"
                for i in range(len(list_of_freq_ptms)+1):
                    ptms_file.write("0\t")
                ptms_file.write("\n")

    print(path_type+" PTMs feature's file created!")


if __name__ == "__main__":

    print("All files are already created and stored in PTMS/ folder \n ")
    while True:
        ans = input("Do you want to recreate files? y/n")
        if ans == "y" or ans == "Y":

            ptmd_file = script_path + "/PTMD.txt"

            list_of_file_names = [fileName for fileName in listdir(script_path + '/data/') if
                                  isfile(join(script_path + "/data/", fileName))]

            list_of_freq_ptms = ["Phosphorylation", "Acetylation", "Ubiquitination", "Succinylation", "Methylation",
                                 "Malonylation", "N-linked Glycosylation", "O-linked Glycosylation", "Sumoylation"]
            #571,032 -> 5,450 experimental sites

            createFeatureFile(list_of_file_names)

            non_pathogenic_feature_file = createAlignmentFile(script_path+"/ptms_feature_file.tsv",
                        script_path+"/filtered_PTMD.tsv", parent_dir_path + "/proteins_alignment/refseq_uniprot.txt")


            createFeatureWithDataset(
                list_of_freq_ptms=list_of_freq_ptms,
                dataset_file= features_dir_path + "/pathogenic_dataset.tsv",
                feature_file= script_path + "/ptms_feature_file.tsv",
                output= "pathogenic_ptms_feature.tsv",
                path_type="pathogenic")

            createFeatureWithDataset(
                list_of_freq_ptms=list_of_freq_ptms,
                dataset_file=features_dir_path + "/non_pathogenic_dataset.tsv",
                feature_file=script_path + "/refseqs_ptms_feature_file.tsv",
                output="non_pathogenic_ptms_feature.tsv",
                path_type= "non pathogenic")

            print("There is no need for merging files automatically since the resulting files are completely"
                  "aligned with the already existing datasets. The merging can be done manually "
                  "by just copying and pasting the *_ptms_features.tsv in the already existing datasets!")

        elif ans == "n" or "N":
            print("Exiting...")
            break
        else:
            print("Wrong inout given, try again by typing y or n.")
