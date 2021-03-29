'''
Necessary files: splitted_data/ folder
'''

import os
from os import listdir
from os.path import isfile, join
import pandas as pd
import re

script_path = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))
parent_dir_path = os.path.abspath(os.path.join(script_path, os.pardir))
features_dir_path = os.path.abspath(os.path.join(parent_dir_path, os.pardir))


def createResultFiles():
	""" Calls isoglypCL.py script for every pathogenic and non pathogenic fasta file and
		produces the resulting files """

	pathogenic_data_path = features_dir_path + "/splitted_data/pathogenic/"
	non_pathogenic_data_path = parent_dir_path + "/splitted_data/non_pathogenic/"

	list_pathogenic_files = [fileName for fileName in listdir(pathogenic_data_path) if isfile(join(pathogenic_data_path, fileName))]
	list_non_pathogenic_files = [fileName for fileName in listdir(non_pathogenic_data_path) if isfile(join(non_pathogenic_data_path, fileName))]
	
	for i in range(1,len(list_pathogenic_files)):
		if list_pathogenic_files[i][-5:] == "fasta":
			print(list_pathogenic_files[i])
			os.system("python3 isoglypCL.py -f " +pathogenic_data_path +"pathogenic_set_final_"+str(i)+".fasta -p isoPara.txt -j /data_files/path_"+str(i))
	
	for i in range(58, len(list_non_pathogenic_files)):
		if list_non_pathogenic_files[i][-5:] == "fasta":
			print(list_non_pathogenic_files[i])
			os.system("python3 isoglypCL.py -f "+ non_pathogenic_data_path+ "non_pathogenic_set_final_"+str(i)+".fasta -p isoPara.txt -j /data_files/non_path_"+str(i))

	print("Results file created!")


def processResults():
	"""
	From each result's file, keep only a value of 0 or 1.
	For each variation we assign:
	-value 1 if the SNP's position is included in the set of the results
	-value 0 else
	"""

	list_results_files = [fileName for fileName in listdir(script_path+"/data_files/") if isfile(join(script_path+"/data_files/", fileName)) and fileName[-4:]==".csv"]
	print(list_results_files)
	for i in range(len(list_results_files)-1, len(list_results_files)):
		if "non_path" in list_results_files[i]: #non pathogenic
			folder_path = script_path +"/non_pathogenic/"
			mut_pos =2

		elif "path" in list_results_files[i]: #pathogenic
			folder_path = script_path + "/pathogenic/"
			mut_pos = 1

		with open(folder_path+list_results_files[i][:-4]+".tsv",'w+') as file_isoglyp: # create file with 1s and 0s
			print(folder_path+list_results_files[i][:-4]+".tsv")
			file_isoglyp.write("Variation"+"\t"+"ISOglyP"+"\n")

			file = pd.read_csv(script_path+"/data_files/"+list_results_files[i]) #, names=["Sequence","S/T","Position","Pattern"]

			prev_variation = file.loc[0][0]
			flag = 0
			for row in range(0, len(file)):
				current_variation = file.loc[row][0]
				mutation = file.loc[row][0].split(" ")[mut_pos]
				position = re.findall(r'\d+', mutation)  # get the positionition of the mutation

				if current_variation != prev_variation:
					write_variation = prev_variation[1:].replace(" ", "|").replace("_",
																				   "|")  # format it like the already existing file
					file_isoglyp.write(write_variation + "\t" + str(flag) + "\n")
					prev_variation = current_variation
					flag =0

				else:
					if int(position[0]) == int(file.loc[row][2]): flag = 1
		break

	print("Results processed successfully")


def mergeFeatureFiles(folder,output,sep=","):
	""" Merges all files of a folder to one file only """

	list_files = [fileName for fileName in listdir(folder) if isfile(join(folder, fileName))]
	df = pd.read_csv(folder+list_files[0], header=None, sep=sep)

	number_of_records = len(df)
	df_complete = df

	for f in range(1, len(list_files)):
		df_new = pd.read_csv(folder+list_files[f], header=None, sep=sep, skiprows=1)
		number_of_records += len(df_new)
		df_complete = pd.concat([df_complete, df_new], ignore_index=True)

	df_complete.to_csv(output, sep='\t', index=False)

	print("Feature merged to a file successfully!")


def mergeFeaturesWithDataset(dataset, feature_file, output, path_type, sep="\t"):
	""" Merges a file which contains a computed feature with one column to the previous complete dataset """

	dataset = pd.read_csv(dataset,sep=sep, low_memory=False, dtype='unicode')
	feature = pd.read_csv(feature_file,sep=sep,skiprows=1)

	labels = dataset["labels"]

	isoglyp_dict = {}
	for label in labels: isoglyp_dict[label] = "None" #dictionary initialization

	if path_type == "non pathogenic":
		for key,value in isoglyp_dict.items(): #fill in with ISOGlyP values
			for feat in feature.values:
				variant = feat[0].split("|")
				protein = variant[2]
				mutation = variant[3][1:] #keep the position and new amino acid to compate

				if key.split("|")[2]==protein and key.split("|")[3][1:] ==mutation :
					isoglyp_dict[key] = feat[1]

	elif path_type == "pathogenic":
		for key,value in isoglyp_dict.items(): #fill in with ISOGlyP values
			for feat in feature.values:
				if key == feat[0]: isoglyp_dict[key] = feat[1]

	dataset["ISOGlyP"] = isoglyp_dict.values()
	out = dataset
	out.to_csv(output, sep="\t",index=False)

	print("Feature merged with "+path_type+" dataset")


def menuMergeChoice():
	"""Providing the choice to the user to merge files since they already exist"""

	while True:
		ans = input("Execute merging with dataset? y/n")

		if ans == "y" or "Y":

			merge_ans = input("merge pathogenic?")
			if merge_ans == "y" or merge_ans == "Y":

				mergeFeaturesWithDataset(dataset=features_dir_path + "/pathogenic_dataset.tsv",
										 feature_file=script_path + "/pathogenic.tsv",
										 output=features_dir_path + "/pathogenic_dataset.tsv",
										 path_type="pathogenic")

			merge_ans = input("merge non pathogenic?")
			if merge_ans == "y" or merge_ans == "Y":

				mergeFeaturesWithDataset(dataset=features_dir_path + "/non_pathogenic_dataset.tsv",
										 feature_file=script_path + "/non_pathogenic.tsv",
										 output=features_dir_path + "/non_pathogenic_dataset.tsv",
										 path_type="non pathogenic")
			break

		elif ans == "n" or "N":
			print("Not merging the feature to file. Exiting script")
			break
		else:
			print("Wrong inout given, try again by typing y or n.")


if __name__ == "__main__":

	print("All files are already created and stored in ISOGlyP-master/ folder \n ")
	while True:
		ans = input("Do you want to recreate files? y/n")
		if ans == "y" or ans == "Y":

			createResultFiles()
			processResults()

			mergeFeatureFiles(script_path + "/pathogenic/", script_path + "/pathogenic.tsv",sep="\t")
			mergeFeatureFiles(script_path + "/non_pathogenic/", script_path + "/non_pathogenic.tsv", sep="\t")

			menuMergeChoice()

		elif ans == "n" or "N":
			print("Will not recreate files! ")
			menuMergeChoice()
			break
		else:
			print("Wrong inout given, try again by typing y or n.")
