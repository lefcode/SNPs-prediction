'''
Necessary files: non_pathogenic_set_final.txt , nps_to_uniprot.tsv and uniprots.txt
'''

import pandas as pd
import os
import csv

script_path = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))
parent_dir_path = os.path.abspath(os.path.join(script_path, os.pardir))
features_dir_path = os.path.abspath(os.path.join(parent_dir_path, os.pardir))
project_dir_path = os.path.abspath(os.path.join(features_dir_path, os.pardir))
dataset_gathering_dir_path = project_dir_path +"/03_Dataset_gathering/data"



def extractProteins(pathogenic_dataset, non_pathogenic_dataset ):
    """Extracts all refseq ids and uniprot ids from the two datasets """

    print("Extracting proteins")

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


def retrieveNPs(non_path_set):
    """Keep only RefSeqs and NPs and store them to a file"""

    print("doing alignment of refseq ids and uniprot ids")

    non_path_file = pd.read_csv(non_path_set, sep=",", names=["rs","mut","np","seq"])

    nps = non_path_file[["rs","np"]]
    nps.to_csv(script_path+"/proteins_ids.txt",index=False,header=None)

    print("alignment done!")
    print("proteins_ids.txt file created!")


def npToUniprot(nps_uni_file, nps_uniprot_file, proteins_ids_file):
    """Align RefSeqs with Uniprot IDs"""

    print("doing alignment with NP ids and uniprot ids")
    
    with open(nps_uni_file,"w+") as nps_uni, open(nps_uniprot_file,"r") as nps_to_uniprot:
        reader = csv.reader(nps_to_uniprot,delimiter='\t' )
        for row in reader:
            nps_uni.write(row[0]+"\t"+row[2]+"\n")

    print("nps_uni.txt file created!")
    rs_nps = pd.read_csv(proteins_ids_file, sep=",", names=["rs","np"])
    np_unis = pd.read_csv(script_path+"/nps_uni.txt" ,sep="\t", names=["nps","uni"])

    print("doing alignment of refseq ids and uniprot ids")
    with open(script_path+"/refseq_uniprot.txt","w") as alignment_file:
        for rs_np in rs_nps.values:
            for np_uni in np_unis.values:
                nps_list = np_uni[0].split(",")
                for i in range(0,len(nps_list)):
                    if nps_list[i] == rs_np[1]:
                        alignment_file.write(rs_np[0] + "," + np_uni[1] + "\n")
    print("alignment done!")
    print("refseq_uniprot.txt file created!")


def uniprotToRefseq(uniprots_file, np_unis_file, rs_nps_file):
    """Align Uniprot IDs with RefSeqs"""

    print("aligning uniprots to refseq ids")

    uniprots = pd.read_csv(uniprots_file)
    np_unis = pd.read_csv(np_unis_file ,sep="\t", names=["nps","uni"])
    rs_nps = pd.read_csv(rs_nps_file, sep=",", names=["rs", "np"])

    with open(script_path+"/uniprot_refseq.txt", 'w+') as alignment_file:
        for uniprot in uniprots.values:
            for np_uni in np_unis.values:
                if uniprot == np_uni[1]: #uniprot found
                    for rs_np in rs_nps.values:
                        if rs_np[1] == np_uni[0]: #same np id
                            # refseq found
                            alignment_file.write(np_uni[1]+","+rs_np[0]+"\n")

    print("uniprot_refseq.txt file created!")
    uniprot_refseq =pd.read_csv(script_path+"/uniprot_refseq.txt", sep=",",names=["uni","ref"])

    pathogenic_refseqs =uniprot_refseq["ref"]
    pathogenic_refseqs.to_csv(script_path+"/pathogenic_refseqs.txt", index=False, header=False)
    print("alignment done!")
    print("pathogenic_refseqs.txt file created")


if __name__ =="__main__":

    extractProteins(features_dir_path + "/pathogenic_dataset/merged_path_dataset.tsv",
                    features_dir_path + "/non_pathogenic_dataset/non_merged_path_dataset.tsv")
                    
    retrieveNPs(dataset_gathering_dir_path+"/non_pathogenic_set_final.txt")
    
    npToUniprot(script_path+"/nps_uni.txt", script_path+"/nps_to_uniprot.tsv",
                script_path+"/proteins_ids.txt")

    uniprotToRefseq(script_path + "/uniprots.txt", script_path + "/nps_uni.txt",
                    script_path + "/proteins_ids.txt")

