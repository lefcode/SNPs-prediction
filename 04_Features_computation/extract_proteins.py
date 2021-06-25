import pandas as pd
import os
import re

script_path = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))

def extractProteinsPositions(pathogenic_dataset, non_pathogenic_dataset ):

    non_path_data = pd.read_csv(non_pathogenic_dataset, sep="\t", low_memory=False)
    labels = non_path_data["labels"]

    with open(script_path+"/refseqs.txt" ,"w+") as refseq_pos:

        for l in labels:
            rs,  = l.split("|")[0]  #keep rs and
            # pos = re.sub("[^0-9]", "", l.split("|")[3])  #position
            refseq_pos.write(rs+ "\n") #","+ str(pos)+


    path_data = pd.read_csv(pathogenic_dataset, sep="\t", low_memory=False)
    labels = path_data["labels"]
    with open(script_path+"/uniprots.txt" ,"w+") as refseq_pos:

        for l in labels:
            uni = l.split("|")[0]
            #pos = re.sub("[^0-9]", "", l.split("|")[1]) # position
            refseq_pos.write(uni+ "\n") #+ str(pos)+


if __name__ == "__main__":

    pathogenic_dataset = script_path + "/new_pathogenic_dataset.tsv"
    non_pathogenic_dataset = script_path + "/new_non_pathogenic_dataset.tsv"

    extractProteinsPositions(pathogenic_dataset, non_pathogenic_dataset)