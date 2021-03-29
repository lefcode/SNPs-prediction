'''
Necessary files: flat chromosomes files, GRCh38_latest_protein.faa and omim files
'''

from os import listdir
from os.path import dirname, isfile, join
import pandas as pd
import gzip
import csv
from Bio import SeqIO

# constants declaration ########################################
VAR = "VAR"
script_path = dirname(__file__)  # gives back your directory path
data_path = script_path + "/data/non_pathogenic_data/"  # path to data
#################################################################


def processFlatFile(chr):
    """
    Implement the necessary conditions on flat files
    We keep only: (missense, validated=YES, notwithdrawn, MAF>0.01, chromosome count>49) records
    From these records we keep the rsID and their locus lines (NP_)
    """

    print("processing chromosome:"+str(chr))
    if chr==4 or chr==5:
        gz = "F:\Διπλωματικη\SNPs Prediction\chromosomes" + "\ds_flat_ch" + str(chr) + ".flat.zip"  # external disk path
    else:
        gz =  "F:\Διπλωματικη\SNPs Prediction\chromosomes"+ "\ds_flat_ch" + str(chr) + ".flat.gz"  # external disk path
    #gz = data_path + "ds_flat_ch" + str(chr) + ".flat.gz" # flat.gz file path
    processed_txt_file = data_path + "filtered_chromosomes/filtered_ch" +str(chr)+".txt"

    with gzip.open(gz, 'r+') as f, open(processed_txt_file,"w+",newline='') as processed_file:
        flag = 0 #flag to check if all conditions are ok
        for line in f:
            if line.strip() and line[:2] != b"ss":  # if line not blank or starting with ss
                row = line.split(b"|")

                if b"rs" in row[0]:  # starts with rsXXXXXX
                    rs_id = str(row[0][:-1].decode()) #keep the rs_id
                    flag =0

                elif b"VAL" in row[0]:  # starts with VAL
                    if (b"validated=YES" in row[1]) and (b"notwithdrawn" in row[4]): # keep if validated and notwithdrawn
                        flag +=1 #ok with this condition

                elif b"GMAF" in row[0]:  # starts with GMAF
                    if (float(row[3][5:]) >= 0.01 and (int(row[2][7:]) > 49)): # keep if MAF>=0.01 and count>49
                        #print(float(row[3][5:8]))
                        flag +=1 #ok with this condition

                elif flag==2 and b'LOC' in row[0]: #all conditions ok
                    try:
                        if (b"fxn-class=missense" in row[3]) and (b'prot_acc=NP_' in row[9]):  # keep if missense and NP_
                            mutation = str(row[4][8:-1].decode() + row[7][13:-1].decode() + row[6][9:-1].decode())  # allele + position + residue
                            processed_file.write(rs_id + ',' + mutation + "," + str(row[9][10:].decode()))
                    except IndexError: #only for chromosome 7, because one record had different format
                        continue #we just ignore that record

    print("chromosome: "+str(chr)+" filtered")


def filterAllFlatFiles():
    """ Calls processFlatFile() to filter all Flat files of all chromosomes"""

    for chr in range(1,23): #do it for all files (1,23)
        processFlatFile(chr)
    processFlatFile("X")
    processFlatFile("Y")


def mergeFilteredFiles(filtered_files_path):
    """ Merges all filtered files to one file complete_filtered.txt"""

    print("merging filtered files")

    #filtered_files_path = data_path + 'filtered_chromosomes'
    list_of_file_paths = [join(filtered_files_path, fileName) for fileName in listdir(filtered_files_path) if
                       isfile(join(filtered_files_path, fileName))]

    df = pd.read_csv(list_of_file_paths[0],header=None,sep=",")
    number_of_records= len(df)
    df_complete =df

    for f in range(1,len(list_of_file_paths)):
        df_new = pd.read_csv(list_of_file_paths[f],header=None,sep=",")
        number_of_records += len(df_new)
        df_complete = pd.concat([df_complete, df_new], ignore_index=True)

    df_complete.to_csv(filtered_files_path +"/complete_filtered.txt" ,sep=',', index=False, header=False)
    print("complete_filtered.txt created!")


def filterOMIM(omim_path, filtered_data):
    """ Filters out all SNPs linked to OMIM """

    print("filtering OMIM data")
    #omim_path = data_path + "omim_files/"
    list_omim_files = [join(omim_path, fileName) for fileName in listdir(omim_path) if
                       isfile(join(omim_path, fileName))]

    new_filtered_file = data_path+"omim_filtered_chromosomes.txt"

    #filtered_data = pd.read_csv(data_path+"/filtered_chromosomes/complete_filtered.txt",sep=",",names=["rsID","mutation","variant"])
    rs_data = filtered_data.drop(["mutation","variant"],axis=1) #keep only rsIDs
    snps_found = 0
    for i,omim_file in enumerate(list_omim_files):
        print("searching in chromosome:"+str(i+1))
        with open(omim_file,'r') as omim_data:
            omim_reader = csv.reader(omim_data,delimiter="\t")
            for i,omim in enumerate(omim_reader):
                if not (len(omim) == 0) and omim[0][0] == "#": #zero elements in the row or not starting with #
                    for rs, index in zip(rs_data.values, rs_data.index.values):
                        if rs[0][2:] == omim[0][1:]:
                            snps_found+=1
                            try:
                                filtered_data.drop(index, inplace=True)
                            except KeyError:  # for the last
                                continue

    print(str(snps_found) + "snps found linked to OMIM ")
    filtered_data.to_csv(new_filtered_file,sep=',', index=False, header=False)


def createFastaSequences(fasta_file):
    """ Creates the fasta sequence for each mutation found in omim_filtered_chromosomes.txt file """

    print("creating fasta sequen")

    #fasta_file = data_path + "GRCh38_latest_protein.faa/GRCh38_latest_protein.faa"
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    mutations_file = data_path+ "omim_filtered_chromosomes.txt"
    mutations_data = pd.read_csv(mutations_file,sep=',',header=None, names=["rsID", "mutation","variant"])

    with open(data_path + "rs_mut_var_fasta.txt", "w+") as rs_mut_var_fasta:
        for f in fasta_dict:
            if not f[:2] == "NP": break #keep only if NP
            for rec in mutations_data.values:
                filter_prot_id = rec[2] # protein in mutations data
                if f == filter_prot_id:  # if found
                    rs_mut_var_fasta.write(rec[0] + "," + rec[1] + "," + rec[2]+ "," +str(fasta_dict[f].seq) + "\n")

    print("fasta sequences created")


def performVariantMutations(variants_file, output_file):
    """ Performs all the mutations to create the fasta sequences of each variant in each protein """

    print("performing the defined mutations in proteins' fasta sequences")
    #variants_file = data_path + "rs_mut_var_fasta.txt"
    #output_file = data_path + "rs_mutations_proteins_fasta_final.txt"
    var_prot_fasta_data = pd.read_csv(variants_file,sep=',',header=None, names=["rsID","mutation", "protein","fasta"])

    with open(output_file,'w+') as rs_variants_proteins_fasta_final:
        for row in var_prot_fasta_data.values:
            after_mut =  row[1][-1] #last character
            mut_position = row[1][1:-1] #4 digits position e.x. 2058

            try: #4 digit position e.x. G3901D
                new_fasta = replaceStrIndex(row[3], int(mut_position)-1, after_mut)
            except ValueError:
                try: #3 digit position e.x. A209S
                    new_fasta = replaceStrIndex(row[3], int(mut_position[0:3])-1, after_mut)
                except ValueError:
                    try: #2 digit position e.x. C64R
                        new_fasta = replaceStrIndex(row[3], int(mut_position[0:2])-1, after_mut)
                    except ValueError: #1 digit position e.x. C5P
                        try:
                            new_fasta = replaceStrIndex(row[3], int(mut_position[0:1])-1, after_mut)
                        except: continue #there was a wrong record which was omitted

            rs_variants_proteins_fasta_final.write(row[0]+","+row[1]+","+row[2]+","+new_fasta+"\n")

    print("mutations applied to fasta sequences")


def replaceStrIndex(text, index, replacement):
    """ Replaces a letter at a position given of a bigger string """

    return '%s%s%s' % (text[:index], replacement, text[index + 1:])


def deleteMultipleRecords(input_file):
    """ Removes records that are identical by comparing the fasta sequences of similar rsIDs that seem to be having
    same mutations """

    #input_file = data_path + "rs_mutations_proteins_fasta_final.txt"
    print("removing duplicates in data")

    output_file = data_path + "non_pathogenic_set.txt"
    file_data = pd.read_csv(input_file,sep=',',header=None, names=["rsID","mutation", "protein","fasta"])

    mutation_fasta_dict ={} #dict
    with open(output_file,'w+') as non_pathogenic_set:
        for rec in file_data.values:
            if not rec[1] in mutation_fasta_dict.keys(): #if not in dictionary,
                mutation_fasta_dict[rec[1]] = rec[3]  #write in dictionary
                new_row =rec.tolist()
                # and in file
                non_pathogenic_set.write(new_row[0]+","+new_row[1]+","+new_row[2]+","+new_row[3]+"\n")

            elif rec[1] in mutation_fasta_dict.keys() : # if mutation in dictionary
                dict_fasta = mutation_fasta_dict.get(rec[1])
                if not (rec[3] ==dict_fasta): #if not the same fasta
                    mutation_fasta_dict[rec[1]] = rec[3]
                    new_row = rec.tolist()
                    non_pathogenic_set.write(new_row[0] + "," + new_row[1] + "," + new_row[2] + "," + new_row[3] + "\n")
                # else: continue # don't write in file

    print("double fasta sequences removed along with their rest records ")


def convert3LetterTo1AminoAcid(input_file):
    """ Converts all mutation 1-letter notations to 3-letter notations so both pathonegic and non-pathogenic files
    have similar format and data """

    print("converting three-letter to one amino acid notation")
    #input_file = data_path + "non_pathogenic_set_final.txt"
    amino_codes_dict = {
        "A": "Αla",  # alanine
        "C": "Cys",  # cystine
        "D": "Asp",  # aspartic acid
        "E": "Glu",  # glutamin acid
        "F": "Phe",  # phenylalanine
        "G": "Gly",  # glycine
        "H": "His",  # histidine
        "I": "Ile",  # isoleucine
        "K": "Lys",  # lysine
        "L": "Leu",  # leucine
        "M": "Met",  # methionine
        "N": "Asn",  # asparagine
        "P": "Pro",  # proline
        "Q": "Gln",  # glutamine
        "R": "Αrg",  # arginine
        "S": "Ser",  # serine
        "T": "Thr",  # threonine
        "V": "Val",  # valine
        "W": "Trp",  # tryptophan
        "Y": "Tyr",  # tyrosine
    }

    data = pd.read_csv(input_file, sep=',', header=None, names=["rsID", "mutation", "protein", "fasta"])
    new_mutations =list()
    for row in data.values:
        mutation = row[1]
        three_letter_before_mutation = amino_codes_dict.get(mutation[0]) #get 3-letter notations
        three_letter_after_mutation =  amino_codes_dict.get(mutation[-1])

        first_replacement = replaceStrIndex(mutation, 0, three_letter_before_mutation)

        new_mutation = replaceStrIndex(first_replacement, len(first_replacement)-1 , three_letter_after_mutation)
        new_mutations.append(new_mutation)

    data['mutation'] = new_mutations
    data.to_csv(input_file,sep=",",index=False, header=False)

    print("1-letter notation replaced by 3-letter notation")


def countProteins(input_file):

    #input_file = data_path+ "non_pathogenic_set_final.txt"
    data = pd.read_csv(input_file,sep=',',header=None,names=["rsID", "mutation", "protein", "fasta"])
    print("mutations:"+ str(len(data)))
    proteins_list = set()
    for row in data.values:
        if row[2][0:6] not in proteins_list:
            proteins_list.add(row[2])

    print("proteins:" + str(len(proteins_list)))


if __name__ =="__main__":

    filterAllFlatFiles()

    mergeFilteredFiles(data_path + 'filtered_chromosomes')

    filterOMIM(data_path + "omim_files/", data_path+"/filtered_chromosomes/complete_filtered.txt")

    createFastaSequences(data_path + "GRCh38_latest_protein.faa/GRCh38_latest_protein.faa")

    performVariantMutations(data_path + "rs_mut_var_fasta.txt", data_path + "rs_mutations_proteins_fasta_final.txt")

    deleteMultipleRecords(data_path + "rs_mutations_proteins_fasta_final.txt")

    convert3LetterTo1AminoAcid( data_path + "non_pathogenic_set_final.txt")

    countProteins(data_path + "non_pathogenic_set_final.txt")
