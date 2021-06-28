# SNPs tools


Part of the total SNPs prediction project <br>
In this part of the project 7 protein features are computed for each pathogenic and non pathogenic SNP. There is no particular order for the features to be computed, except from the DEOGEN2 which requires InterProScan folder's data. <br>
In order to use the repository locally follow the next steps

## Open folder and execute the corresponding script for the tool you want to use

`cd SNPs_tools/`

#### proteins_alignment folder is already complete. refseq_to_uniprot.py script does every necessary alignment between uniprot ids and refseq ids. It can be executed using the following command:

`python3 proteins_alignment/refseq_to_uniprot.py`


### 1. ISOGlyP
This tool predicts if a spesific position of a protein is glycosylated and creates a single-column feature with 1 and 0:

`python3 InterProScan/isoglyp.py`

### 2. InterProScan
This folder examines the domains that a protein contains and creates an equivalent file for each dataset with columns for the 10 most frequent and one for the non frequent domains:

`python3 InterProScan/interpro.py`

### 3. DEOGEN2
This part computes the PFscore and creates a single-column feature with numerical continuous values:

`python3 DEOGEN2/deogen.py`

### 4. MVP
This folder examines whether a protein takes part in complex formations based on CORUM database and creates a single-column feature with 1 and 0:

`python3 MVP/mvp.py`

### 5. PTMs: 
Based on PTMD and dbPTM databases which contain information for ptms happening in spesific positions of proteins, we find the 9 most f1equent ptms and create a 0/1 feature for each one of them and one feature representing all non frequent ptms:

`python3 PTMs/ptms.py`

### 6. Q(SASA) and structure annotations

`python3 ZoomVar_query_script/zoomvar.py`
