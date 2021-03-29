'''
Necessary files: CORUM_allComplexes.txt
http://mips.helmholtz-muenchen.de/corum/#download
'''

import pandas as pd
import os

script_path = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))
dataPath = script_path + "/data/"
parent_dir_path = os.path.abspath(os.path.join(script_path, os.pardir))
features_dir_path = os.path.abspath(os.path.join(parent_dir_path, os.pardir))


def complexFormationIdentification(pathogenic_proteins_file, non_pathogenic_proteins_file, corum_file):
    """Since CORUM database includes proteins that are present in complex formations,
    we find out, for each of our proteins, if they are included in this file.
    We firstly filter out any proteins that are found in other species rather than humans and create a
    file assigning 1 or 0 depending on protein existance or not in corum file"""

    print("identifying proteins included in CORUM dataset. (complex formations)")

    corum_data = pd.read_csv(corum_file, sep="\t").drop(
        ["ComplexID", "ComplexName", "Cell line", "subunits(Entrez IDs)",
         "Protein complex purification method", "GO ID", "GO description",
         "FunCat description", "subunits(Gene name)", "Subunits comment",
         "Complex comment", "SWISSPROT organism", "subunits(Gene name syn)",
         "subunits(Protein name)", "PubMed ID", "Disease comment",
         "FunCat ID"], axis=1)

    corum_data.drop(corum_data[corum_data['Organism'] != "Human"].index, inplace=True)

    pathogenic_proteins = pd.read_csv(pathogenic_proteins_file)
    non_pathogenic_proteins = pd.read_csv(non_pathogenic_proteins_file, sep=',',
                                          names=["rs", "uni"])

    with open(script_path + "/path_corum.txt", "w+") as path_file:
        path_file.write("Uniprot\tComplexFormation\n")
        for row in corum_data.values:
            proteins_list = row[2].split(";")
            for prot in proteins_list:
                if prot in pathogenic_proteins.values:  
                    path_file.write(prot + "\t" + str(1)+"\n")
                else:
                    path_file.write(prot + "\t" + str(0)+"\n")

    print("path_corum.txt file created!")

    with open(script_path + "/non_path_corum.txt", "w+") as non_path_file:
        non_path_file.write("RefSeq\tComplexFormation\n")

        for rs_uni in non_pathogenic_proteins.values:
            complex_flag = 0
            for row in corum_data.values:
                proteins_list = row[2].split(";")

                if rs_uni[1] in proteins_list:
                    complex_flag = 1
                    non_path_file.write(rs_uni[0] + "\t" + str(1) + "\n")
                    break

            if complex_flag == 0:
                non_path_file.write(rs_uni[0] + "\t" + str(0) + "\n")

    print("non_path_corum.txt file created!")


def mergeFeaturesWithDataset(dataset, feature_file, output, path_type, feature_name, sep="\t"):
    """ Merges a file which contains a computed feature with one column to the previous complete dataset"""

    dataset = pd.read_csv(dataset, sep=sep, low_memory=False, dtype='unicode')
    feature = pd.read_csv(feature_file, sep=sep, skiprows=1)

    labels = dataset["labels"]

    data_dict = {}
    for label in labels: data_dict[label] = "0"  # dictionary initialization

    if path_type == "non pathogenic":
        for key in data_dict.keys():
            for feat in feature.values:
                protein = feat[0]
                if protein == key.split("|")[0]:  # when find the same protein id then asign the computed values
                    data_dict[key] = feat[1]

    elif path_type == "pathogenic":
        for key, value in data_dict.items():  # fill in with ISOGlyP values
            for feat in feature.values:
                protein = feat[0]
                if protein == key.split("|")[0]:
                    data_dict[key] = feat[1]

    dataset[feature_name] = data_dict.values()
    out = dataset
    out.to_csv(output, sep="\t", index=False)

    print("Feature merged with " + path_type + " dataset")


def menuMergeChoice():
    """Providing the choice to the user to merge files since they already exist"""

    while True:
        ans = input("Execute merging with dataset? y/n")
        if ans == "y" or "Y":

            merge_ans = input("merge pathogenic?")
            if merge_ans == "y" or merge_ans == "Y":

                mergeFeaturesWithDataset(dataset=features_dir_path + "/pathogenic_dataset.tsv",
                                         feature_file=script_path + "/path_corum.txt",
                                         output=features_dir_path + "/pathogenic_dataset.tsv",
                                         feature_name="ComplexFormation", path_type="pathogenic")

            merge_ans = input("merge non pathogenic?")
            if merge_ans == "y" or merge_ans == "Y":

                mergeFeaturesWithDataset(dataset=features_dir_path + "/non_pathogenic_dataset.tsv",
                                         feature_file=script_path + "/non_path_corum.txt",
                                         output=features_dir_path + "/non_pathogenic_dataset.tsv",
                                         feature_name="ComplexFormation", path_type="non pathogenic")
            break

        elif ans == "n" or "N":
            print("Not merging the feature to file. Exiting script")
            break
        else:
            print("Wrong inout given, try again by typing y or n.")


if __name__ == "__main__":

    print("All files are already created and stored in MVP/data/ and MVP/ folders \n ")
    while True:
        ans = input("Do you want to recreate files? y/n")
        if ans == "y" or ans == "Y":

            corum_file = dataPath + "CORUM_allComplexes.txt"

            pathogenic_proteins_file = parent_dir_path + "/proteins_alignment/uniprots.txt"
            non_pathogenic_proteins_file = parent_dir_path + "/proteins_alignment/refseq_uniprot.txt"

            complexFormationIdentification(pathogenic_proteins_file, non_pathogenic_proteins_file, corum_file)

            menuMergeChoice()

        elif ans == "n" or "N":
            print("Will not recreate files! ")
            menuMergeChoice()
            break
        else:
            print("Wrong inout given, try again by typing y or n.")