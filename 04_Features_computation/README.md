# Features Computation


This part of the project receives the two types of data: pathogenic and non pathogenic protein sequences in fasta format and explores their possible features.
In tools folder new features are computed based on computational tools and public databases' filtering.

The new data are merged to the previously computed features that lie into pathogenic_dataset and non_pathogenic_dataset folders. These features were computed in [EnsembleGASVR]
(https://www.researchgate.net/publication/261922806_EnsembleGASVR_A_novel_ensemble_method_for_classifying_missense_Single_Nucleotide_Polymorphisms)
and are recomputed and used in this project



1. O-Glycosylation

This feature refers to the prediction whether the position of a variant is glycosylated. The computation of this feature involves downloading the [ISOGlyP tool] (https://github.com/jonmohl/ISOGlyP) which predicts if a serine or threonine is glycosylated in all positions of the protein. ISOGlyP is followed by [ISOGlyP-EV_Tables] (https://github.com/jonmohl/ISOGlyP-EV_Tables) that should be downloaded and stored at the same directory as ISOGlyP. Fistly, we devide the data to smaller datablocks and execute the following command for all fasta files.
```
python3 isoglypCL.py -f /path/to/file.fasta -p isoPara.txt -j output_file
```
From the resulting csv files we check for each variant whether the position of the mutation is included in the "Position" column so we receive a 1 for the variant otherwise a 0.



2. Protein Domains

In this feature we examine whether each variant includes one (or many) of the 10 most frequent protein domains or any of the rest non-frequent. To compute this feature we downloaded from [InterProScan](https://www.ebi.ac.uk/interpro/download/) the interpro_entries_list.tsv
and from [UniProtKB](https://www.uniprot.org/uniprot/?query=&sort=score) the uniprot_interpro.tsv file by maintaining only Entry and Cross-reference (InterPro) columns. Since in the non pathogenic dataset we have refseq ids we firstly had to do an alignment between Uniprot ids and refseq ids using [Retrieve/ID mapping tool of UniProt] (https://www.uniprot.org/uploadlists/)

The procedure is as follows. We align the uniprot and refseq ids with IPR identifiers of InterPro domains and replace the IPRs with the domain names. Then we count the 10 most freqeunt domain names found in the alignment files. Lastly, we create 10 columns with the 10 most frequent domains and 1 which corresponds to all the rest domains together. For each protein we assign 1 to the corresponding column if the protein included each domain, else 0. For the non frequent domains feature, if the protein includes any of the non frequent domains then we assign 1, else 0.



3. PF score

This feature is closely related to the previous feature and quantifies the sensitivity of domains to pathogenic variations.
For each domain we count how many neutral and pathogenic variants correspond to them and compute the PF score based on the following formula:

pf_score_domain = log((Ndel +1)/(Ndel + Nneut +2)) - log((Nneut +1)/(Ndel + Nneut + 2))

where Ndel = # of deleterious variants and Nneut = # of neutral variants

Then, for each variant we assign the PF score of the domain that it includes.
If a variant includes multiple domains, we assign the maximum value amongst the PF scores.
If a variant doesn't include a domain then we assign 0 value. 



4. Complex Formation

In this feature we find all proteins that are involved into protein complexes. We downloaded from [CORUM database](http://mips.helmholtz-muenchen.de/corum/#download) the CORUM_allComplexes.txt file and we filtered this file to maintain only records from the "subunits(UniProt IDs)" column that refer to Human. For each protein of our datasets we checked whether it is involved in any protein complex and we assigned 1 to the ones that do, otherwise we assigned 0.


5. PTMs

In this feature we find whether a variant falls into a same position with frequent and non-frequent PTMs.
We downloaded from [PTMD database](http://ptmd.biocuckoo.org/download.php) the PTMD.txt file and from [dbPTM database](http://dbptm.mbc.nctu.edu.tw/download.php) all gzipped files.
We firstly filter out all the non-human records from the PTMD file and each of the dbPTM files and maintain the uniprot id, the position and the name of the PTM. Afterwards we find the 9 most frequent PTMs .
Then, we make an alignment between uniprot and refseq ids to create a pathogenic and a non-pathogenic files with the PTMs and 
Lastly, for each variant of the dataset we inspect whether the position falls into a frequent PTM, so we assign 1 to the equivalent column, else 0 or falls into any of the non-frequent PTMs and the 10th column gets 1, else 0



6. ZoomVar features ( Qsasa and structure annotations)

This feature computes the Q(SASA) and structure annotations based on ZoomVar tool. Before the execution an alignment must take place because the ZoomVar tool requires refseq ids as input. Thus, we created a refseq id list that corresponds to the pathogenic data.

In order for this folder to work appropriately the following steps must be followed:

1. Install [ZoomVar Rest API](http://fraternalilab.kcl.ac.uk/ZoomVar/downloads/)
2. Install [VEP command line tool](https://www.ensembl.org/info/docs/tools/vep/script/index.html)
3. Install full cache file [homo_sapiens_vep_102_GRCh38.tar.gz](ftp://ftp.ensembl.org/pub/current_variation/indexed_vep_cache/)
4. Execute VEP commands for non pathogenic and pathogenic ids:
```
./vep --cache --format id --symbol --canonical --protein -i ~/link/to/SNPs-Prediction/04_Features_computation/proteins_alignment/refseqs.txt -o non_pathogenic_vep.txt --force

./vep --cache --format id --symbol --canonical --protein -i ~/link/to/SNPs-Prediction/04_Features_computation/proteins_alignment/pathogenic_refseqs.txt -o pathogenic_vep.txt --force
```

In case of error with refseq.txt and pathogenic_refseqs.txt files: move them into VEP_input_data and execute again the previous commands by giving the new path of these files as input.

5. Execute ZoomVar Rest API commands for structure information only:
```
python3 rest_script.py ~/link/to/SNPs_tools/ZoomVar_query_script/ensembl-vep/non_pathogenic_vep.txt non_pathogenic_structure.csv -v -s -n 20 -i 20
python3 rest_script.py ~/link/to/SNPs_tools/ZoomVar_query_script/ensembl-vep/pathogenic_vep.txt pathogenic_structure.csv -v -s -n 20 -i 20
```

Subsequently, for each variant we assign the resulting values that ZoomVar produces and are stored to non_pathogenic_structure.csv and pathogenic_structure.csv files. The assignment happens by checking the position and the name of the variant.
Since the resulting files have multiple records for some variants, but with different values we assign the maximum of these values.

