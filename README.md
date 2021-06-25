# SNPs prediction Integrated Msc thesis

The purpose of this work is to predict whether a given SNP (Single Nucleotid Polymorphism) is deleterious or neutral.

The complete work is devided in 3 parts.


# Dataset gathering

In order to develop the project and create a classifier to predict SNP's pathogenity both pathogenic and benign data were collected.


Pathogenic data were retrieved from the online database [PhenCode](http://phencode.bx.psu.edu/) and non pathogenic data from
[dbSNP](https://www.ncbi.nlm.nih.gov/snp/) <br>

Regarding PhenCode database the appropriate filtering was then done in order to keep only non-synomynous SNPs related to deseases affecting the phenotype. <br>

dbSNP data were also processed to keep only missense SNPs with MAF>0.01, count>49 and notwithdrawn existing to real proteins. Proteins connected to [OMIM database](https://www.omim.org/) were also removed from our collection.  <br>

For both datasets we needed to retrieve the fasta sequences for each protein so we used the [Retrieve/ID](https://www.uniprot.org/uploadlists/) tool of [UniProt](https://www.uniprot.org/) and the [RefSeq FTP](https://ftp.ncbi.nlm.nih.gov/refseq/) tool of [NCBI](https://www.ncbi.nlm.nih.gov/)

Subsequently, we filter all the common sequences that exist in both sets, because they cannot be labeled simultaneously as pathogenic and neutral.

Finally, we collect with similar procedure the data from HumVar and HumDiv datasets from [Polyphen-2] 
(http://genetics.bwh.harvard.edu/pph2/dokuwiki/downloads) in order to use them for comparisons.


# Features computation

The purpose of this step was to retrieve information i.e. features from the proteins that we collected. Many features were already computed and mined from T.Rapakoulia in [EnsemblGASVR](https://www.researchgate.net/publication/261922806_EnsembleGASVR_A_novel_ensemble_method_for_classifying_missense_Single_Nucleotide_Polymorphisms) so the main goal was for this dataset to be integrated and enriched.

Information was found for O-glycosylation, Protein Domains and their sensitivity towards pathogenity due to mutations, PTMs, Solvent Accessibility Surface Area and Protein-Protein Interactions.


Finally, an 18 features dataset 105 columns) was produced for training and testing.


# Training Classifiers


As a last step, the data from the previous step are preproccessed and given as input to an optimized Random Forest Classifier, while the most important features are extracted.

The ability of the classifier to make correct predictions is evaluated with a separate and distinct dataset, which was extracted from the collected data, as well as with the HumVar and HumDiv datasets

The results are very satisfying and promising for SNPs prediction

