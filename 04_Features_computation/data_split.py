'''
Necessary files:  pathogenic_set_final.fasta , non_pathogenic_set_final.fasta , pathogenic_set_final.txt
and non_pathogenic_set_final.txt
'''

import pandas as pd
from Bio import SeqIO
from os.path import dirname,abspath,join
import os
from fsplit.filesplit import Filesplit
fs = Filesplit()


script_path = dirname(__file__)  # gives back your directory path
splitted_data_path = script_path+"/splitted_data/"
parent_path = abspath(join(script_path, os.pardir))
data_path = parent_path + "/03_Dataset_gathering/data/"  # path to data


def splitPathogenicData(pathogenic_file, length):
    fs.split(file=pathogenic_file, split_size=length, output_dir=splitted_data_path+"pathogenic/")


def splitNonPathogenicData(non_pathogenic_file, length):
    fs.split(file=non_pathogenic_file, split_size=length, output_dir=splitted_data_path+"non_pathogenic/")


def findSimilarMutations(pathogenic_file_txt,non_pathogenic_file_txt):
    """ Finds all similar mutations to the dataset created by John Thanos and merges the data """

    prev_pathogenic = script_path+"/pathogenic_dataset/merged_path_dataset.tsv"
    prev_non_pathogenic = script_path + "/non_pathogenic_dataset/merged_non_path_dataset.tsv"

    prev_pathogenic_data = pd.read_csv(prev_pathogenic,sep="\t",dtype=str)
    prev_non_pathogenic_data = pd.read_csv(prev_non_pathogenic, sep="\t",dtype=str)

    pathogenic_data = pd.read_csv(pathogenic_file_txt,sep=',',header=None,
                                  names=["mutation", "protein","fasta"]).drop(["mutation","fasta"],axis=1)
    non_pathogenic_data = pd.read_csv(non_pathogenic_file_txt,sep=',',header=None,
                                    names=["rsID", "mutation", "protein", "fasta"]).drop(["protein","fasta"],axis=1)


    with open(script_path+"/new_pathogenic_set.tsv",'w+') as new_pathogenic_set:
        new_pathogenic_set.write('\t'.join((list(map(str, prev_pathogenic_data.columns.tolist())))) + "\n")
        for row in pathogenic_data.values:
            for prev in prev_pathogenic_data.values:
                if row[0][:6] == prev[0][:6] and row[0][10:]==prev[0][8:]:
                    #print('\t'.join((list(map(str, prev.tolist()))))+"\n")
                    new_pathogenic_set.write('\t'.join((list(map(str, prev.tolist()))))+"\n")
    
    print("pathogenic set finished")

    with open(script_path+"new_non_pathogenic_set.tsv", 'w+') as new_non_pathogenic_set:
        for row in non_pathogenic_data.values:
            for prev in prev_non_pathogenic_data.values:
                print(prev[0].split("|"))
                #if row[0][:6] == prev[0][:6]:
                #    new_non_pathogenic_set.write('\t'.join((list(map(str, prev.tolist()))))+"\n")

    print("non pathogenic set finished")


def extractLabels(pathogenic_file_fasta,non_pathogenic_file_fasta):
    pathogenic_dict = SeqIO.to_dict(SeqIO.parse(pathogenic_file_fasta, "fasta"))
    non_pathogenic_dict = SeqIO.to_dict(SeqIO.parse(non_pathogenic_file_fasta, "fasta"))


if __name__=="__main__":
    pathogenic_file_fasta = data_path + "pathogenic_set_final.fasta"
    non_pathogenic_file_fasta = data_path + "non_pathogenic_set_final.fasta"

    pathogenic_file_txt = data_path+"pathogenic_set_final.txt"
    non_pathogenic_file_txt = data_path+"non_pathogenic_set_final.txt"

    #extractLabels(pathogenic_file_fasta,non_pathogenic_file_fasta)
    #findSimilarMutations(pathogenic_file_txt,non_pathogenic_file_txt)

    length = 5000000  # 900000 -> 1200 records
    splitPathogenicData(pathogenic_file_fasta,length)
    splitNonPathogenicData(non_pathogenic_file_fasta,length)