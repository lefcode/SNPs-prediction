# SNPs_tools

Part of the total SNPs prediction project. In order to use the repository locally follow the next steps

## Clone repository
`git clone https://github.com/lefcode/SNPs_tools.git`

## Open folder and execute the corresponding script for the tool you want to use

`cd SNPs_tools/`

#### proteins_alignment folder is already complete. refseq_to_uniprot.py script does every necessary alignment between uniprot ids and refseq ids. It can be executed using the following command:

`python3 proteins_alignment/refseq_to_uniprot.py`


### 1. for DEOGEN2, part which computes the pfscore and creates a single-column feature with numerical continuous values:

`python3 DEOGEN2/deogen.py`

### 2. for InterProScan, folder which examines the domains that a protein contains and creates an equivalent file for each dataset with columns for the 10 most frequent and one for the non frequent domains:

`python3 InterProScan/interpro.py`

### 3. for ISOGlyP, tool which predicts if a spesific position of a protein is glycosylated and creates a single-column feature with 1 and 0:

`python3 InterProScan/isoglyp.py`

### 4. for MVP, folder which examines whether a protein takes part in complex formations based on CORUM database and creates a single-column feature with 1 and 0:

`python3 MVP/mvp.py`

### 5. for PTMs: based on PTMD and dbPTM databases. which contain information for ptms happening in spesific positions of proteins, we find the 9 most frequent ptms and create a 0/1 feature for each one of them and one feature representing all non frequent ptms:

`python3 PTMs/ptms.py`

### 6. for ZoomVar, tool which computes the Q(SASA) and the number of PPIs happening in spesific positions of proteins.

In order for this script to run it is necessary for VEP command line tool
( https://www.ensembl.org/info/docs/tools/vep/script/index.html ) to be installed along with full cache file: 
homo_sapiens_vep_102_GRCh38.tar.gz (ftp://ftp.ensembl.org/pub/current_variation/indexed_vep_cache/)

For the pathogenic file an alignment between uniprot ids and refseq ids was needed in order to run the VEP command line tool since it requires only refseq ids as input. This is executed already in proteins_alignment/refseq_to_uniprot.py

ZoomVar script requires a VEP file so both files must be given to rest_api.py as input.

### VEP commands
```
./vep --cache --format id --symbol --canonical --protein -i ~/link/to/SNPs_Tools/ZoomVar_query_script/VEP_input_data/refseqs.txt -o non_pathogenic_vep.txt --force
./vep --cache --format id --symbol --canonical --protein -i ~/link/to/SNPs_Tools/ZoomVar_query_script/VEP_input_data/pathogenic_refseqs.txt -o pathogenic_vep.txt --force
```

### After executing all the previous the REST API can now be executed. use the following commands:
```
python3 rest_script.py ~/link/to/SNPs_tools/ZoomVar_query_script/ensembl-vep/non_pathogenic_vep.txt non_pathogenic_structure.csv -v -s -n 20 -i 20
python3 rest_script.py ~/link/to/SNPs_tools/ZoomVar_query_script/ensembl-vep/pathogenic_vep.txt pathogenic_structure.csv -v -s -n 20 -i 20
```

### Lastly, the following command will construct the feature columns and merge them to the initial datasets
`python3 ZoomVar_query_script/zoomvar.py`
