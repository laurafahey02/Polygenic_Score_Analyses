<!-- dx-header -->
#  PRSice-2

## What does this app do?

PRSice (pronounced 'precise') is a Polygenic Risk Score app for calculating, applying, evaluating and plotting the results of polygenic risk scores (PRS) analyses.

## What are typical use cases for this app?

The app is used for calculating, applying, evaluating, and plotting the results of polygenic risk scores (PRS) analyses.
PRSice2 can run at high resolution to provide the best-fit PRS as well as provide results calculated at broad P-value thresholds, illustrating results corresponding to either, can thin SNPs according to linkage disequilibrium and P-value ("clumping"), and can be applied across multiple traits in a single run.

## What data are required for this app to run?

PRSice-2 accepts:

- A base data file `*.assoc` containing GWAS summary statistics. This file should contain the following columns (please note that the exact header names given below are expected by the app):
    - SNP: SNP ID, usually in the form of rs-ID.
    - CHR: The chromosome in which the SNP resides.
    - BP: Chromosomal co-ordinate of the SNP.
    - A1: The effect allele of the SNP.
    - A2: The non-effect allele of the SNP.
    - P: The P-value of association between the SNP genotypes and the base phenotype.
    - OR: The effect size estimate of the SNP, if the outcome is binary/case-control. If the outcome is continuous or treated as continuous then this will be the BETA.
- A PLINK BED file `*.bed` containing the primary representation of genotype calls at biallelic variants in the BED format. Please note that the file should only contain biallelic variants.
- A PLINK BIM file `*.bim` is the extended variant information file accompanying a .bed binary genotype table. Please note that the GWAS summary base file coordinates and target geno files (bed,bim) should be of same genome build.
- A PLINK FAM file `*.fam` contains the sample information file accompanying a .bed binary genotype table.
- (Optional) A  phenotype `.pheno` file containing the phenotype, covariates of the samples including the principal components (PC). Please note that this input is beneficial to test more than one phenotype and when including covariates. The file should contain the following columns:
    - FID: Family identifier
    - IID: Individual identifier
    - Phenotype (e.g. height)
- A trait choice on whether the trait is 'binary' value or 'quantitative'.
- (Optional) A boolean plot generation option that generates a bar plot with PRS results corresponding to a range of P-value thresholds and a scatter plot with "best-fit" PRS and the phenotype.

## What does this app output?

A `.PRSice` text file containing all the polygenic risk scores, a `.summary` summary text file containing the summary about the phenotype, SNPs, polygenic scores, and other metrics, and a `.best` text file containing the best polygenic scores that are not overfitting to the model based on this tool’s prediction. Optionally, two additional files can be returned: one containing a bar plot with P-value threshold vs PRS model fit and one with a scatter plot with PRS vs Phenotype.

## How does this app work?

This app runs a Docker image containing the PRSice-2 package. It utilizes the base file provided to calculate PRS and generate metrics to understand the model fit.

Documentation can be found at https://www.prsice.info/
