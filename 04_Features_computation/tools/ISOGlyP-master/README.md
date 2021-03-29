# ISOGlyP
Isoform Specific mucin-type o-glycosylation prediction

# Parameter File
Modify the parameter file to update THR/SER ratios, transferases desired and positions

# Cite ISOGlyP
Coming soon!

# Usage
```
isoglypCL.py [-h] [-f FASTA] [-a ACCESSIONS] [-p PARAMETERS] [-j JOBID]
                    [-c CUTOFF] [-l WRITE_LOG]

optional arguments:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
  -a ACCESSIONS, --accessions ACCESSIONS
  -p PARAMETERS, --parameters PARAMETERS
  -j JOBID, --jobId JOBID
  -c CUTOFF, --cutoff CUTOFF
  -l WRITE_LOG, --write_log WRITE_LOG
```


Run for our data

python3 isoglypCL.py -f ~/Desktop/SNPs_Prediction/03_Dataset_gathering/data/test.fasta -p isoPara.txt -j path1

python3 isoglypCL.py -f ~/Desktop/SNPs_Prediction/03_Dataset_gathering/data/pathogenic_set_final.fasta -p isoPara.txt
python3 isoglypCL.py -f ~/Desktop/SNPs_Prediction/03_Dataset_gathering/data/non_pathogenic_set_final.fasta -p isoPara.txt

