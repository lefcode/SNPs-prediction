## ZoomVar features ( Q(SASA) and structure interactions) computation

In order for this folder to work appropriately the following steps must be followed:

#### 1. Install [ZoomVar Rest API](http://fraternalilab.kcl.ac.uk/ZoomVar/downloads/)
#### 2. Install [VEP command line tool](https://www.ensembl.org/info/docs/tools/vep/script/index.html)
#### 3. Install full cache file [homo_sapiens_vep_102_GRCh38.tar.gz](ftp://ftp.ensembl.org/pub/current_variation/indexed_vep_cache/)
#### 4. Execute VEP commands:
```
./vep --cache --format id --symbol --canonical --protein -i ~/link/to/SNPs-Prediction/04_Features_computation/proteins_alignment/refseqs.txt -o non_pathogenic_vep.txt --force
./vep --cache --format id --symbol --canonical --protein -i ~/link/to/SNPs-Prediction/04_Features_computation/proteins_alignment/pathogenic_refseqs.txt -o pathogenic_vep.txt --force
```

In case of error with refseq.txt and pathogenic_refseqs.txt files: move them into VEP_input_data 
and execute again the previous commands by giving the new path of these files as input.

#### 5. Execute ZoomVar Rest API commands:
```
python3 rest_script.py ~/link/to/SNPs_tools/ZoomVar_query_script/ensembl-vep/non_pathogenic_vep.txt non_pathogenic_structure.csv -v -s -n 20 -i 20
python3 rest_script.py ~/link/to/SNPs_tools/ZoomVar_query_script/ensembl-vep/pathogenic_vep.txt pathogenic_structure.csv -v -s -n 20 -i 20
```

### 6. Execute the following command to construct the feature columns and merge them to the datasets
```
python3 zoomvar.py
```