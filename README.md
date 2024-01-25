# Polygenic Score Analyses
Analyses performed for the study, "Shared Genetic Links between Sleep, Neurodevelopmental, and Neuropsychiatric Conditions: A Pathway-Based Polygenic Score Analysis".


## Preparation of Target Data

### UK Biobank QC on the UKB RAP

#### Step 1: Sample QC
UKB Research participants were first excluded based on the following criteria: not being of white British ancestry (data-field: 22006), containing chromosomal aneuploidies (data-field: 22019), having a high SNP missingness and/or having unusually high or low heterozygosity (data-field: 22027) and performing shift work (data-field: 826). These samples were identified by UKB (Bycroft et al., 2018) and accessible via the data fields reported. 

Participants are furthermore removed based on relatedness using the UKB supplied relatedness file, which lists pairs of individuals related up to a third degree. One individual from each pair was removed without removing samples from pairs where a sample had already been removed.

This was done using Python in jupyter notebooks on UKB RAP as detailed here: create_list_samples_to_remove.ipynb

#### Step 2: Create list SNPs pass the imputation information threshold
Create a list of SNPs to be excluded based on imputation quality score (INFO) < 0.7.

This was done using R in jupyter notebooks on UKB RAP as detailed here: CreateListSNPsPassImpINFO.ipynb

#### Step 3: Run Plink on genotype files to exclude SNPs and samples that do not pass quality checks
Bash command

## Preparation of discovery data
### For genome-wide analyses using [SBayesRC](https://www.biorxiv.org/content/10.1101/2022.10.12.510418v1)


### For pathway-based analyses using [PRSet](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010624) 
