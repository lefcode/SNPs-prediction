from os.path import dirname
import pandas as pd
from non_pathogenic_preprocess import convert3LetterTo1AminoAcid

script_path = dirname(__file__)  # gives back your directory path
data_path = script_path + "/data/"  # path to data


def removeCommonSequences(pathogenic_file,non_pathogenic_file,output_pathogenic_file,output_non_pathogenic_file):
    """ Finds and removes all fasta sequences and their equivalent record that exist in both
    pathogenic and non pathogenic datasets """

    print("removing common sequences of two datasets")

    pathogenic_data= pd.read_csv(pathogenic_file,sep=',',header=None, names=["mutation", "protein","fasta"])
    non_pathogenic_data = pd.read_csv(non_pathogenic_file,sep=',',header=None,names=["rsID","mutation", "protein","fasta"])
    count_equals = 0 #count the equal fasta sequences

    for row in non_pathogenic_data.values:
        ind = pathogenic_data.index[pathogenic_data['fasta'] == row[3]].tolist() # get all indices where fasta sequences are equal
        if not len(ind)==0: #if found any
            count_equals+= len(ind)
            pathogenic_data.drop(ind,inplace=True) #drop
            non_pathogenic_data.drop(ind, inplace=True) #drop

    print(str(count_equals) + "sequences found and removed from both datasets")
    pathogenic_data.to_csv(output_pathogenic_file,sep=',', index=False, header=False)
    non_pathogenic_data.to_csv(output_non_pathogenic_file, sep=',', index=False, header=False)


def convertToFasta(input_file, choice):
    """ Converts a txt file to .fasta file
    :param input_file: The file
    :param choice: choice between pathogenic (0) and non pathogenic (1)
    """

    print("converting file to .fasta ")

    if choice ==0:
        input_data = pd.read_csv(input_file,sep=',',header=None,names=["mutation", "protein","fasta"])
        with open(input_file[:-4] + ".fasta", 'w+') as file:
            for i in range(0, len(input_data)):
                protein = input_data.loc[i]["protein"]
                file.write(">"+ protein[:6] +" "+ protein[9:] +"\n" + input_data.loc[i]["fasta"] +"\n" )

    elif choice == 1:
        input_data = pd.read_csv(input_file, sep=',', header=None, names=["rsID","mutation", "protein", "fasta"])
        with open(input_file[:-4] + ".fasta", 'w+', encoding="utf-8") as file:
            for i in range(0, len(input_data)):
                file.write(">"+ input_data.loc[i]["rsID"]+ " " + input_data.loc[i]["mutation"] + " " +
                           input_data.loc[i]["protein"] +"\n"+input_data.loc[i]["fasta"] +"\n" )
    else:
        print('wrong input parameters')
        return 0

    print('file converted to fasta file')


if __name__ == "__main__":
    pathogenic_file= data_path + "pathogenic_data/pathogenic_set.txt"
    output_pathogenic_file = data_path +"pathogenic_data/pathogenic_set_final.txt"

    non_pathogenic_file= data_path + "non_pathogenic_data/non_pathogenic_set.txt"
    output_non_pathogenic_file =  data_path + "non_pathogenic_data/non_pathogenic_set_final.txt"

    removeCommonSequences(pathogenic_file,non_pathogenic_file,output_pathogenic_file,output_non_pathogenic_file)

    choice = 0
    convertToFasta(data_path + "pathogenic_data/pathogenic_set_final.txt", choice)

    choice = 1
    convertToFasta(data_path + "non_pathogenic_data/non_pathogenic_set_final.txt", choice)
