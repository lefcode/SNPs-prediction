'''
Necessary files: gv.txt , gvAttr.csv and uniprot.sprot.fasta
'''

from os.path import dirname
import pandas as pd
from Bio import SeqIO

# constants declaration ##########################################
VAR = "VAR"
script_path = dirname(__file__)  # gives back your directory path
data_path = script_path + "/data/pathogenic_data/"
##################################################################

def readGVData(gv_file, gv_attr_file):
    """ Stores gv and gvAttr data as pandas dataframes and calls processGVData() """

    gv_data = pd.read_csv(gv_file, sep="\t", header=None,
                         names=["mutation", "protein", "polymorphism", "type", "sequence", "number"])
    gv_attr_data = pd.read_csv(gv_attr_file, sep=",",header=None,encoding='utf-8')
    gv_attr_data.columns=["mutation","process","effect"]

    return gv_data,gv_attr_data


def writeFilteredDataToTxt(gv_data, gv_attr_data, output_file):
    """ Produces the filtered_data.txt which contains 19.861 variants and 1.743 proteins,
        based on the mutation that appear in both gv and gvAttr files """

    print("writing filtered data to filtered_data.txt ")
    #output_file =data_path+"filtered_data.txt"
    with open(output_file,'w+',newline='') as filtered_file:
        #csv_writer = csv.writer(file, delimiter=',') #, quotechar='"', quoting=csv.QUOTE_MINIMAL
        for row in gv_data.values:
            for sp in gv_attr_data.values:
                if row[0] == sp[0]:
                    filtered_file.write(row[0]+","+row[1]+"\n")


def processGVData(gv_file, gv_attr_file):
    """ Filters the gv and gvAttr files based on our conditions to create the filtered_data.txt """

    print("processing gv files")
    gv_data,gv_attr_data = readGVData(gv_file, gv_attr_file)

    """ For the gv.txt we keep only, SPs, substitutions and changes in 1 amino acid"""
    print("Initial gv.txt size:", len(gv_data))
    gv_data.drop(gv_data[gv_data['polymorphism'] != "SP"].index, inplace=True)
    gv_data.drop(gv_data[gv_data['type'] != "substitution"].index, inplace=True)
    gv_data.drop(gv_data[gv_data['number'] != 1].index, inplace=True)
    gv_data.drop(gv_data[gv_data['protein'].str[:6] == "Q9Y2S0"].index, inplace=True) #proteins to be removed
    print("Size of gv.txt after pre-processing:", len(gv_data))

    """ For the gvAttr.txt we keep only phenotype-associated diseases and only Variants"""
    print("Initial gvAttr.csv size:", len(gv_attr_data))
    gv_attr_data.drop(gv_attr_data[gv_attr_data['effect'] != "phenotype-associated"].index, inplace=True)
    gv_attr_data.drop(gv_attr_data[gv_attr_data['process'] != "disease"].index, inplace=True)
    gv_attr_data.drop(gv_attr_data[gv_attr_data['mutation'].str[:3] != VAR ].index, inplace=True)
    print("Size of gvAttr.csv after pre-processing:", len(gv_attr_data))

    writeFilteredDataToTxt(gv_data,gv_attr_data)


def maintainAllIdentifiers(filtered_file, output_file):
    """ Keeps only the protein identifiers in order to create the fasta sequences.
     Creates identifiers.csv file which later will be used to find the fasta sequences of each protein """

    print("maintaining all proteins identifiers (uniprot ids) and stores them to identifiers.txt")
    # filtered_file = data_path+"filtered_data.txt"
    # output_file = data_path+'identifiers.txt'
    filtered_data = pd.read_csv(filtered_file, sep=",", header=None,
                         names=["mutation", "protein"])

    identifiers = filtered_data['protein']
    identifiers_filtered = list()
    for id in identifiers.values:
        identifiers_filtered.append(id[0:6])

    identifiers_filtered_frame = pd.DataFrame(identifiers_filtered)

    identifiers_filtered_frame.to_csv(output_file,index=False)


def createVariantsFastas(uniprot_fasta_file, filtered_file, output_file):
    """ Creates the var_prot_fasta.txt file which contains all the variants, proteins: mutation and fasta
    sequences of each protein. calls performVariantMutations() which will contain all the """

    print("creating var_prot_fasta.txt file ")

    #uniprot_fasta_file = data_path + "/uniprot_sprot.fasta/uniprot_sprot.fasta"
    # filtered_file = data_path + "filtered_data.txt"
    # output_file = data_path+ "var_prot_fasta.txt"
    uniprot_fasta_dict = SeqIO.to_dict(SeqIO.parse(uniprot_fasta_file, "fasta"))

    filtered = pd.read_csv(filtered_file, sep=',',header=None, names=["mutation", "protein"])

    with open(output_file ,"w+") as var_prot_fasta:
        for f in uniprot_fasta_dict:
            fasta_prot_id = f[3:9] # protein identifier

            for rec in filtered.values:
                filter_prot_id = rec[1][0:6] # protein in filtered data
                if fasta_prot_id == filter_prot_id: # if found
                    var_prot_fasta.write(rec[0]+ "," + rec[1] +","+ str(uniprot_fasta_dict[f].seq)+"\n")


def performVariantMutations(variants_file, output_file):
    """ Performs all the mutations to create the fasta sequences of each variant in each protein """

    print("performing the defined mutations to the fasta sequences")
    #variants_file = data_path + "var_prot_fasta.txt"
    #output_file =data_path+"pathogenic_set.txt"
    var_prot_fasta_data = pd.read_csv(variants_file,sep=',',header=None, names=["mutation", "protein","fasta"])

    amino_codes_dict = {
        "A": "Αla", #alanine
        "B": "Asp", #aspartate,asparagine
        "C": "Cys", #cystine
        "D": "Asp", #aspartic acid
        "E": "Glu", #glutamin acid
        "F": "Phe", #phenylalanine
        "G": "Gly", #glycine
        "H": "His", #histidine
        "I": "Ile", #isoleucine
        "K": "Lys", #lysine
        "L": "Leu", #leucine
        "M": "Met", #methionine
        "N": "Asn", #asparagine
        "P": "Pro", #proline
        "Q": "Gln", #glutamine
        "R": "Αrg", #arginine
        "S": "Ser", #serine
        "T": "Thr", #threonine
        "V": "Val", #valine
        "W": "Trp", #tryptophan
        "Y": "Tyr", #tyrosine
    }

    with open(output_file,'w+') as variants_proteins_fasta_final:
        for row in var_prot_fasta_data.values:
            #before_mut = row[1][9:12]
            after_mut =  row[1][-3:]
            mut_position = row[1][12:16] #4 digits position e.x. 2058

            for key, value in amino_codes_dict.items():
                if after_mut == value:
                    amino_after = key

            try: #4 digit position e.x. P35498:p.Leu1309Phe
                new_fasta = replaceStrIndex(row[2], int(mut_position)-1, amino_after)
            except ValueError:
                try: #3 digit position e.x. P35498:p.Glu954Lys
                    new_fasta = replaceStrIndex(row[2], int(mut_position[0:3])-1, amino_after)
                except ValueError:
                    try: #2 digit position e.x. P81172:p.Cys70Arg
                        new_fasta = replaceStrIndex(row[2], int(mut_position[0:2])-1, amino_after)
                    except ValueError: #1 digit position e.x. 7
                        new_fasta = replaceStrIndex(row[2], int(mut_position[0:1])-1, amino_after)

            variants_proteins_fasta_final.write(row[0]+","+row[1]+","+new_fasta+"\n")


def replaceStrIndex(text,index,replacement):
    """ Replaces a letter at a position given of a bigger string """

    return '%s%s%s'%(text[:index],replacement,text[index+1:])


def countProteins(input_file):
    """Counts all distinct proteins of the pathogenic dataset """

    #input_file = data_path+ "non_pathogenic_set_final.txt"
    data = pd.read_csv(input_file,sep=',',header=None,names=["mutation", "protein", "fasta"])
    print("mutations:"+ str(len(data)))
    proteins_list = set()
    for row in data.values:
        if row[1][0:6] not in proteins_list:
            proteins_list.add(row[1][0:6])

    print("proteins:" + str(len(proteins_list)))


if __name__ =="__main__":
    print("processing starts!")

    processGVData(data_path + "gv.txt", data_path + "gvAttr.csv")

    maintainAllIdentifiers(data_path+"filtered_data.txt", data_path+'identifiers.txt')

    createVariantsFastas(data_path + "/uniprot_sprot.fasta/uniprot_sprot.fasta",
                         data_path + "filtered_data.txt", data_path+ "var_prot_fasta.txt")

    performVariantMutations(data_path + "var_prot_fasta.txt",data_path + "pathogenic_set.txt")

    countProteins(data_path+"pathogenic_set_final.txt")