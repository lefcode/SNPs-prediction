'''
Necessary files: InterPro files
'''

import os
import pandas as pd
from math import log

script_path = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))
parent_dir_path = os.path.abspath(os.path.join(script_path, os.pardir))
features_dir_path = os.path.abspath(os.path.join(parent_dir_path, os.pardir))


def countDomains(path_domains_file, non_path_domains_file, domains_file):
    """ Since InterProScan folder contains pathogenic_domains.tsv and non_pathogenic_domains.tsv
    we use the available information in order to create domains_count.tsv file. For each domain we count
    all the pathogenic and non pathogenic proteins they it's included into. """

    print("counting neutral and deleterious proteins per domain")

    path_domains = pd.read_csv(path_domains_file, sep="\t")
    non_path_domains = pd.read_csv(non_path_domains_file, sep="\t")
    domains = pd.read_csv(domains_file, sep="\t", names =["Domain","counts"]).drop(["counts"],axis=1)

    with open(script_path+"/pfscore/domains_count.tsv", "w+") as dom_file:
        dom_file.write("Domain name\tneutral;deleterious\n")

        for d in domains.values:
            domain = d[0]
            domain_path_count = 0
            domain_non_path_count = 0

            for uni_dom in non_path_domains.values:
                if type(uni_dom[2]) != float: # if not nan
                    domain_list = uni_dom[2].split(";")
                    for dom in domain_list:
                        if dom == domain:
                            domain_non_path_count += 1

            for uni_dom in path_domains.values:
                if type(uni_dom[1]) != float:  # if not nan
                    domain_list = uni_dom[1].split(";")

                    for dom in domain_list:
                        if dom == domain:
                            domain_path_count+=1

            dom_file.write(domain+"\t"+str(domain_path_count)+";"+str(domain_non_path_count)+"\n")

    print("domains_count.tsv file created")


def domainPFscore(domains_count):
    """ Computes the pf scores for each domain (Domain log-odd score)
    We compute the PF score for each domain which is derived from Deogen2 supplementary material.
    The score is the following:
    PFscore(d) = log((Ndel +1)/(Ndel + Nneut +2)) - log((Nneut +1)/(Ndel + Nneut + 2))
    """

    print("computing pf score of each domain")
    domains_count = pd.read_csv(domains_count, sep="\t")
    with open(script_path+"/pfscore/domains_pfscores.tsv","w+") as domains_pfscores:
        for d in domains_count.values:

            Nneut = int(d[1].split(";")[0])
            Ndel = int(d[1].split(";")[1])

            pf_score_d = log((Ndel +1)/(Ndel + Nneut +2)) - log((Nneut +1)/(Ndel + Nneut + 2))

            domains_pfscores.write(d[0]+"\t"+str(pf_score_d)+"\n")

    print("domains_pfscores.tsv file created!")


def createFeatureFiles(domains_pfscores, path_domains_file, non_path_domains_file):
    """ Creating feature file with intention to merge it to the dateset"""

    print("creating feature files with pf score for each protein")
    path_domains = pd.read_csv(path_domains_file, sep="\t")
    non_path_domains = pd.read_csv(non_path_domains_file, sep="\t")
    pfscores = pd.read_csv(domains_pfscores, sep="\t")

    with open(script_path+"/pfscore/pathogenic_pfscore_feature.tsv","w+") as path_feature_file:
        path_feature_file.write("Uniprot\tPFscore\n")
        for path in path_domains.values:

            if type(path[1]) != float: #not nan
                domains_list = path[1].split(";")

                pfs_list = list()
                for dom in domains_list:
                    for pfscore in pfscores.values:
                        if pfscore[0] == dom:
                            pfs_list.append(float(pfscore[1]))
                            
                try:
                    path_feature_file.write(path[0]+"\t"+str(max(pfs_list))+"\n")
                except ValueError:
                    path_feature_file.write(path[0]+"\t"+"0"+"\n")

            else:
                path_feature_file.write(path[0]+"\t"+"0"+"\n")

    print("pathogenic feature file created!")

    with open(script_path+"/pfscore/non_pathogenic_pfscore_feature.tsv","w+") as non_path_feature_file:
        non_path_feature_file.write("RefSeq\tPFscore\n")
        for non_path in non_path_domains.values:

            if type(non_path[2]) != float: #not nan
                domains_list = non_path[2].split(";")

                pfs_list = list()
                for dom in domains_list:
                    for pfscore in pfscores.values:
                        if pfscore[0] == dom:
                            pfs_list.append(float(pfscore[1]))

                try:
                    non_path_feature_file.write(non_path[0]+"\t"+str(max(pfs_list))+"\n")
                except ValueError:
                    non_path_feature_file.write(non_path[0]+"\t"+"0"+"\n")
            else:
                non_path_feature_file.write(non_path[0]+"\t"+"0"+"\n")

    print("non pathogenic feature file created!")


def mergeFeatureWithDataset(dataset_file, feature_file, output, path_type, sep="\t"):
    """Merges a file which contains the computed features, i.e. 1 column to the previous complete dataset"""

    print("Merging feature with " + path_type + " dataset")

    dataset = pd.read_csv(dataset_file, sep=sep, low_memory=False)
    feature = pd.read_csv(feature_file, sep=sep)

    labels = dataset["labels"]

    complex_dict = {}  # values will be replaced with the computed values, else they remain None if not found
    for label in labels: complex_dict[label] = "None"  # dictionary initialization

    # fill in with complex formation values
    for key in complex_dict.keys():
        for feat in feature.values:
            protein = feat[0]
            if protein == key.split("|")[0]:  # when find the same protein id then asign the computed values
                complex_dict[key] = feat[1]

    dataset["PFscore"] = complex_dict.values()
    out = dataset
    out.to_csv(output, sep="\t", index=False)

    print("PF score feature merged to " + path_type + "dataset!")


def menuMergeChoice():
    """Providing the choice to the user to merge files since they already exist"""

    while True:
        ans = input("Execute merging with dataset? y/n")
        if ans == "y" or ans == "Y":

            merge_ans = input("merge pathogenic?")
            if merge_ans == "y" or merge_ans == "Y":

                mergeFeatureWithDataset(
                    dataset_file=features_dir_path + "/pathogenic_dataset.tsv",
                    feature_file=script_path + "/pfscore/pathogenic_pfscore_feature.tsv",
                    output=features_dir_path + "/pathogenic_dataset.tsv",
                    path_type="pathogenic")

            merge_ans = input("merge non pathogenic?")
            if merge_ans == "y" or merge_ans == "Y":

                mergeFeatureWithDataset(
                    dataset_file=features_dir_path + "/non_pathogenic_dataset.tsv",
                    feature_file=script_path + "/pfscore/non_pathogenic_pfscore_feature.tsv",
                    output=features_dir_path + "/non_pathogenic_dataset.tsv",
                    path_type="non pathogenic")
            break

        elif ans == "n" or "N":
            print("Not merging the feature to file. Exiting script")
            break
        else:
            print("Wrong inout given, try again by typing y or n.")


if __name__ =="__main__":

    print("All files are already created and stored in pfscore/ folder \n ")
    while True:
        ans = input("Do you want to recreate files? y/n")
        if ans == "y" or ans == "Y":
            domains_count_file = script_path + "/pfscore/"

            countDomains(parent_dir_path+"/InterProScan/data/pathogenic_domains.tsv",
                         parent_dir_path+"/InterProScan/data/non_pathogenic_domains.tsv",
                         parent_dir_path+"/InterProScan/domains.txt")

            domainPFscore(script_path+"/pfscore/domains_count.tsv")

            createFeatureFiles(script_path+"/pfscore/domains_pfscores.tsv",
                               parent_dir_path+"/InterProScan/data/pathogenic_domains.tsv",
                               parent_dir_path+"/InterProScan/data/non_pathogenic_domains.tsv")

            menuMergeChoice()

        elif ans == "n" or "N":
            print("Will not recreate files! ")
            menuMergeChoice()
            break
        else:
            print("Wrong inout given, try again by typing y or n.")
