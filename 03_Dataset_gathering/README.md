# Data set gathering

This part of the project is essential in order to collect all the necessary data for the analysis.
SNPs can be categorized into to categories with respect to their pathogenity. A SNP can be either pathogenic or non-pathogenic.  
Thus, it's vital for us to gather data from both categories of SNPs.

The files used in this part of the project can be found and downloaded from the following links.

## Pathogenic data

PhenCode database is used for the pathogenic data.
gv.txt and gvAttr.txt files were downloaded from http://phencode.bx.psu.edu/dist/phencode/database/
using the following commands:
```
wget http://phencode.bx.psu.edu/dist/phencode/database/gv.txt
wget http://phencode.bx.psu.edu/dist/phencode/database/gvAttr.txt 
```
Subsequently, all fasta sequences of the proteins were needed for the analysis.
UniProt database provides a fasta file with reviewed Swiss-Prot proteins and it's available 
[here](https://www.uniprot.org/downloads).

This file was also downloaded to produce the final pathogenic dataset. 

In order to create the pathogenic dataset the user must execute the following command:
```
python3 pathogenic_preprocess.py
```

## Non Pathogenic data

dbSNP database is used for the non-pathogenic data.
The flat files of all chromosomes (1-22, X and Y) were downloaded from 
[dbSNP FTP Download section](https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/ASN1_flat/).

One can use the following command:
```
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/ASN1_flat/ds_flat_ch<chromosome#>.flat.gz
```
Apart from the flat files we need the fasta sequences of the proteins, so we downloaded them from
the [RefSeq FTP server](ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_protein.faa.gz)

One can use the following command:
```
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_protein.faa.gz 
```

Lastly, we downloaded the OMIM files from internet in order to remove any proteins related to the OMIM database.

All mentioned files were downloaded and used from the equivalent script.
Their sizes are prohibitive to be uploaded in this repository, thus any user wanting to 
run the script must first download them and put them in the same folder.

In order to create the pathogenic dataset the user must execute the following command:
```
python3 non_pathogenic_preprocess.py
```

Here, we remind that this is not needed and the final non-pathogenic dataset is provided in the folder.


## Comparison between data sets

As a last step, all common sequences must be removed from both datasets and files should be converted to fasta files in 
order to be used in the next step of features computation.

The equivalent script can be executed using the following command 
```
python3 datasets_comparison.py
```
