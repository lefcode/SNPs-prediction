'''
Necessary files: interpro_entries_list.tsv and uniprot_interpro.tsv
'''

import pandas as pd
import os

script_path = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))
parent_dir_path = os.path.abspath(os.path.join(script_path, os.pardir))
features_dir_path = os.path.abspath(os.path.join(parent_dir_path, os.pardir))
project_dir_path = os.path.abspath(os.path.join(features_dir_path, os.pardir))
dataset_gathering_dir_path = project_dir_path +"/03_Dataset_gathering/data"


def extractProteins(protein_column):
    """Extract all uniprot ids and store them to a file"""

    print("extracting uniprot ids")
    uniprots = [prot[0][:6] for prot in protein_column.values]
    with open(parent_dir_path + "/proteins_alignment/uniprots.txt", 'w') as file:
        for uni in uniprots:
            file.write('%s\n' % uni)
    return uniprots


def filterUniIPR(uni_inter_file, pathogenic_file, non_pathogenic_file, output_pathogenic_file,
                 output_non_pathogenic_file):
    """ Aligns all uniprot ids with the IPR (interpro) ids and stores them to 2 different files
    depending on mutation's pathogenity or neutrality """

    uni_inter = pd.read_csv(uni_inter_file, sep="\t", skiprows=1, names=["Uniprot", "InterPro"])

    pathogenic = pd.read_csv(pathogenic_file, sep=",", header=None, names=["var", "prot", "seq"]).drop(["var", "seq"],
                                                                                                       axis=1)
    non_pathogenic = pd.read_csv(non_pathogenic_file, sep=",", header=None, names=["rs", "uni"])

    uniprots = extractProteins(pathogenic)
   
    with open(output_pathogenic_file, "w+") as path_file:
        path_file.write("Uniprot\tIPR\n")
        for uni in uniprots:
            for row in uni_inter.values:
                if uni == row[0]:
                    path_file.write(uni + "\t" + str(row[1]) + "\n")
    print("Pathogenic uniprot ids are now aligned with IPRS!")
    print(output_pathogenic_file+" file created!")

    with open(output_non_pathogenic_file, "w+") as non_path_file:
        non_path_file.write("RefSeq\tUniprot\tIPR\n")
        for ref_uni in non_pathogenic.values:
            for row in uni_inter.values:
                if ref_uni[1] == row[0]:
                    non_path_file.write(ref_uni[0] + "\t" + ref_uni[1] + "\t" + str(row[1]) + "\n")

    print("Non pathogenic refseqs and uniprot ids are now aligned with IPRS!")
    print(output_non_pathogenic_file+" file created!")


def alignDomains(interpro_entries_list, pathogenic_uni_ipr, non_pathogenic_uni_ipr,
                 output_pathogenic_file, output_non_pathogenic_file):
    """ Creates 2 files containing all uniprot ids, iprs and the domain names. Again 2 files depending on
    mutation's pathogenity or neutrality """

    print("Aligning domains to create pathogenic and non pathogenic feature files")

    interpro_entries = pd.read_csv(interpro_entries_list, sep="\t")
    interpro_entries = interpro_entries.drop(interpro_entries[interpro_entries["ENTRY_TYPE"] != "Domain"].index)

    pathogenic = pd.read_csv(pathogenic_uni_ipr, sep="\t")  # maybe drop_duplicates() needed
    non_pathogenic = pd.read_csv(non_pathogenic_uni_ipr, sep="\t")

    with open(output_pathogenic_file, "w+") as path_file:
        path_file.write("Uniprot\tDomains\n")
        for path in pathogenic.values:
            path_file.write(path[0] + "\t")
            if type(path[1]) != float:  # if not nan
                for inter in interpro_entries.values:
                    ipr_list = path[1].split(";")  # all iprs
                    for ipr in ipr_list:
                        if inter[0] == ipr:
                            path_file.write(inter[2] + ";")

                path_file.write("\n")

            else:
                path_file.write("None\n")
    
    print("Pathogenic file ready! Domains are now aligned with uniprot ids and IPRs")
    print(output_pathogenic_file + " file created!")

    with open(output_non_pathogenic_file, "w+") as non_path_file:
        non_path_file.write("RefSeq\tUniprot\tDomains\n")
        for non_path in non_pathogenic.values:
            non_path_file.write(non_path[0] + "\t" + non_path[1] + "\t")
            if type(non_path[2]) != float:  # if not nan
                ipr_list = non_path[2].split(";")  # all iprs
                for inter in interpro_entries.values:
                    for ipr in ipr_list:
                        if inter[0] == ipr:
                            non_path_file.write(inter[2] + ";")

                non_path_file.write("\n")

            else:
                non_path_file.write("None\n")

    print("Non pathogenic file ready! Domains are now aligned with refseqs, uniprot ids and IPRs")
    print(output_non_pathogenic_file + " file created!")


def findMostFrequentDomains(pathogenic_domains_file, non_pathogenic_domains_file):
    """ Finds the 10 most frequent domains in both pathogneic and non_pathogenic files """

    print("finding 10 most frequent domains")
    pathogenic_domains = pd.read_csv(pathogenic_domains_file, sep="\t")
    non_pathogenic_domains = pd.read_csv(non_pathogenic_domains_file, sep="\t")

    # pathogenic domains
    path_domains_dict = {}
    for path in pathogenic_domains.values:
        if type(path[1]) != float:
            domains_list = path[1].split(";")
            for dom in domains_list:
                if (dom != ""):
                    if (dom not in path_domains_dict):
                        path_domains_dict[dom] = 1
                    else:
                        path_domains_dict[dom] += 1

    non_path_domains_dict = {}
    for non_path in non_pathogenic_domains.values:
        if type(non_path[2]) != float:
            domains_list = non_path[2].split(";")
            for dom in domains_list:
                if (dom != ""):
                    if (dom not in non_path_domains_dict):
                        non_path_domains_dict[dom] = 1
                    else:
                        non_path_domains_dict[dom] += 1

    # create a signle dictionary containing all domains
    for key, value in non_path_domains_dict.items():
        if key not in path_domains_dict:
            path_domains_dict[key] = value
        else:
            path_domains_dict[key] += value  # if domain in both then add the two values

    path_domains_list = sorted(path_domains_dict.items(), key=lambda item: item[1], reverse=True)
    with open(script_path + "/domains.txt", "w+") as domains_file:
        domains_file.write('\n'.join('{}\t{}'.format(x[0], x[1]) for x in path_domains_list))  # write tuple

    print("domains.txt file created!")

    domains = list() # keep only the domains and not their appearance numbers
    for row in path_domains_list:
        domains.append(row[0])

    domains.remove("Immunoglobulin-like domain")
    domains.remove("Immunoglobulin subtype")

    frequent_domains_list = [domains[i] for i in range(10)]
    print("10 most frequent domains found!")
    print(frequent_domains_list)

    non_frequent_domains_list = [domains[i] for i in range(10,len(domains))]

    return frequent_domains_list, non_frequent_domains_list


def createDomainFeatures(frequent_domains_list, non_frequent_domains_list, pathogenic_domains, non_pathogenic_domains,
                         output_pathogenic_feats, output_non_pathogenic_feats):
    """ Asserts 1 if the domain is included in the most frequent otherwise 0 """

    pathogenic_domains = pd.read_csv(pathogenic_domains, sep="\t")
    non_pathogenic_domains = pd.read_csv(non_pathogenic_domains, sep="\t")

    print("creating domains feature files")

    with open(output_pathogenic_feats, "w+") as path_file:
        path_file.write("Uniprot\t")
        for i in range(len(frequent_domains_list)):
            path_file.write(frequent_domains_list[i] + "\t")
            # path_file.write("Domain{}\t".format(str(i))) #writing the title
        path_file.write("Non-Freq-Domain\n")  # +1 due to the extra non frequent domain feature

        for path in pathogenic_domains.values:
            path_file.write(path[0] + "\t")
            domain_flag = 0
            if type(path[1]) != float:
                dom_list = path[1].split(";")
                for i in range(10):  # iterate 11 times for all domains + the non frequent domain
                    # print(path_frequent_domains[i])
                    if (dom_list != "") and (frequent_domains_list[i] in dom_list):
                        path_file.write("1\t")
                    else:
                        path_file.write("0\t")

                non_frequent_domain_flag = 0
                for d in dom_list:
                    if d in non_frequent_domains_list:

                        non_frequent_domain_flag =1

                if non_frequent_domain_flag == 0:
                    path_file.write("0\n")
                elif non_frequent_domain_flag == 1:
                    path_file.write("1\n")

            else:
                for i in range(10):
                    path_file.write("None\t")

                path_file.write("None\n")

            '''
            if domain_flag == 0:  # Non frequent domain
                path_file.write("1\n")
            elif domain_flag == 1:
                path_file.write("0\n")
            else:
                path_file.write("None\n")
            '''

    print(output_pathogenic_feats + " file created!")

    with open(output_non_pathogenic_feats, "w+") as non_path_file:
        non_path_file.write("RefSeq\tUniprot\t")
        for i in range(len(frequent_domains_list)):
            non_path_file.write(frequent_domains_list[i] + "\t")  # writing the title
        non_path_file.write("Non-Freq-Domain\n")  # +1 due to the extra non frequent domain feature

        for non_path in non_pathogenic_domains.values:
            non_path_file.write(non_path[0] + "\t" + non_path[1] + "\t")
            domain_flag = 0
            if type(non_path[2]) != float:
                dom_list = non_path[2].split(";")
                for i in range(10):  # iterate 11 times for all domains + the non frequent domain
                    # print(path_frequent_domains[i])
                    if (dom_list != "") and (frequent_domains_list[i] in dom_list):
                        non_path_file.write("1\t")
                    else:
                        non_path_file.write("0\t")

                non_frequent_domain_flag = 0
                for d in dom_list:
                    if d in non_frequent_domains_list:
                        non_frequent_domain_flag = 1

                if non_frequent_domain_flag == 0:
                    non_path_file.write("0\n")
                elif non_frequent_domain_flag == 1:
                    non_path_file.write("1\n")

            else:
                for i in range(10):
                    non_path_file.write("None\t")

                non_path_file.write("None\n")

            '''
            if domain_flag == 0:  # Non frequent domain
                non_path_file.write("1\n")
            elif domain_flag == 1:
                non_path_file.write("0\n")
            else:
                non_path_file.write("None\n")
            '''
    print(output_non_pathogenic_feats + " file created!")

    print("Domain features created!")
    return frequent_domains_list


def mergeFeaturesWithDataset(domain_list, dataset_file, feature_file, output, path_type, sep="\t"):
    """Merges a file which contains the computed features, i.e. 11 columns to the previous complete dataset"""

    dataset = pd.read_csv(dataset_file, sep=sep, low_memory=False)
    features = pd.read_csv(feature_file, sep=sep)

    labels = dataset["labels"]

    if path_type == "non pathogenic":
        protein_id = "RefSeq"  # indicating the id to be compared, RefSeq of Uniprot
    elif path_type == "pathogenic":
        protein_id = "Uniprot"

    print("merging feature files with " + path_type + " dataset!")

    for i in range(11):  # iterate 11 times for all domains +the non frequent domain
        domain_dict = {}  # values will be replaced with the computed values, else they remain None if not found
        for label in labels: domain_dict[label] = "None"  # dictionary initialization

        if i == 10: domain_name = "Non-Freq-Domain"  # last domain
        else: domain_name = domain_list[i]

        print(domain_name)
        # fill in with Domain values
        for key in domain_dict.keys():
            for feat in features[[protein_id, domain_name]].values:
                # only the specific domain-column
                protein = feat[0]
                if protein == key.split("|")[0]:  # when find the same protein id then asign the computed values
                    domain_dict[key] = feat[1]

        dataset[domain_name] = domain_dict.values()
        out = dataset
        out.to_csv(output, sep="\t", index=False)

    print(path_type + " dataset merge with features")


def menuMergeChoice(frequent_domains_list):
    """Providing the choice to the user to merge files since they already exist"""

    while True:
        ans = input("Execute merging with dataset? y/n")
        if ans == "y" or "Y":

            merge_ans = input("merge pathogenic?")
            if merge_ans == "y" or merge_ans == "Y":

                mergeFeaturesWithDataset(
                    domain_list=frequent_domains_list,
                    dataset_file=features_dir_path + "/pathogenic_dataset.tsv",
                    feature_file=script_path + "/pathogenic_feats.tsv",
                    output=features_dir_path + "/pathogenic_dataset.tsv",
                    path_type="pathogenic")

            merge_ans = input("merge non pathogenic?")
            if merge_ans == "y" or merge_ans == "Y":

                mergeFeaturesWithDataset(
                    domain_list=frequent_domains_list,
                    dataset_file=features_dir_path + "/non_pathogenic_dataset.tsv",
                    feature_file=script_path + "/non_pathogenic_feats.tsv",
                    output=features_dir_path + "/non_pathogenic_dataset.tsv",
                    path_type="non pathogenic")
            break

        elif ans == "n" or "N":
            print("Not merging the feature to file. Exiting script")
            break
        else:
            print("Wrong inout given, try again by typing y or n.")


if __name__ == "__main__":

    print("All files are already created and stored in InterProScan/data/ and InterProScan/ folders \n ")
    while True:
        ans = input("Do you want to recreate files? y/n")

        if ans == "y" or ans == "Y":
            '''
            filterUniIPR(script_path + "/data/uniprot_interpro.tsv",
                        dataset_gathering_dir_path +"/pathogenic_set_final.txt",
                        parent_dir_path + "/proteins_alignment/refseq_uniprot.txt",
                        script_path + "/data/pathogenic_uni_ipr.tsv",
                        script_path + "/data/non_pathogenic_uni_ipr.tsv")

            interpro_entries_list = script_path + "/data/interpro_entries_list.tsv"
            
            alignDomains(interpro_entries_list,
                         script_path + "/data/pathogenic_uni_ipr.tsv",
                         script_path + "/data/non_pathogenic_uni_ipr.tsv",
                         script_path + "/data/pathogenic_domains.tsv",
                         script_path + "/data/non_pathogenic_domains.tsv")
            '''

            frequent_domain_list, non_frequent_domains_list = findMostFrequentDomains(
                                                    script_path + "/data/pathogenic_domains.tsv",
                                                    script_path + "/data/non_pathogenic_domains.tsv")

            createDomainFeatures(frequent_domain_list, non_frequent_domains_list,
                                 script_path + "/data/pathogenic_domains.tsv",
                                 script_path + "/data/non_pathogenic_domains.tsv",
                                 script_path+"/pathogenic_feats.tsv",
                                 script_path+"/non_pathogenic_feats.tsv")

            menuMergeChoice(frequent_domain_list)

        elif ans == "n" or "N":
            print("Will not recreate files! ")
            frequent_domain_list, non_frequent_domains_list = findMostFrequentDomains(
                                    script_path + "/data/pathogenic_domains.tsv",
                                    script_path + "/data/non_pathogenic_domains.tsv")
            menuMergeChoice(frequent_domain_list)
            break
        else:
            print("Wrong inout given, try again by typing y or n.")
