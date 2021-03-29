# SNPs prediction Integrated Msc thesis

The purpose of this work is to predict whether a given SNP (Single Nucleotid Polymorphism) is deleterious or neutral.

The complete work is devided in 5 parts.

# Biological background
Firstly, there is a chapter explaining the fundamentals that this thesis is based on i.e. DNA, mutation's basics, SNP definition etc.

# Existing methodologies research

Secondly, a folder containing many published papers-tools is created. An extended research was done to gain insight into other techniques and methods including various types of data being analyzed as well as for comperative purposes.

# Dataset gathering

As a third part of this project many data were retrieved from the online databases [PhenCode](http://phencode.bx.psu.edu/) and [dbSNP](https://www.ncbi.nlm.nih.gov/snp/). 
Regarding PhenCode database the appropriate filtering was then done in order to keep only non-synomynous SNPs related to deseases affecting the phenotype. 
dbSNP data were also processed to keep only missense SNPs with MAF>0.01, count>49 and notwithdrawn existing to real proteins. Proteins connected to [OMIM database](https://www.omim.org/) were also removed from our collection.
For both datasets we needed to retrieve the fasta sequences for each protein so we used the [Retrieve/ID](https://www.uniprot.org/uploadlists/) tool of [UniProt](https://www.uniprot.org/) and the [RefSeq FTP](https://ftp.ncbi.nlm.nih.gov/refseq/) tool of [NCBI](https://www.ncbi.nlm.nih.gov/). This step was vital for the next step.

# Features computation

The purpose of this step was to retrieve information i.e. features from the proteins that we collected. Many features were already computed and mined from T.Rapakoulia [EnsemblGASVR](https://www.researchgate.net/publication/261922806_EnsembleGASVR_A_novel_ensemble_method_for_classifying_missense_Single_Nucleotide_Polymorphisms) so the main goal was for this dataset to be integrated and enriched.

[ISOGlyP](https://isoglyp.utep.edu/) tool was used to predict if a protein is glycosylated or not.

Additionaly, we retrieved the 10 most frequent protein domains found in [InterProScan](https://www.ebi.ac.uk/interpro/search/sequence/). An extra feature is dedicated to non-frequent domains.

Whether a protein has a complex formation is another vital information. Thus, we used [CORUM database](http://mips.helmholtz-muenchen.de/corum/) which provides this information.


Another feature proposed in [DEOGEN2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5570203/) is PFscore, a score to quantify the domain information that a protein lies into.


