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
The following code was run from the UKB RAP command line interface.
```bash
for chr in {1..22}; do
   plink_qc="plink2 --bgen project-GP8V3yjJXgZZ9vbKBFz74V2Q:/Bulk/Imputation/UKB imputation from genotype/ukb22828_c${chr}_b0_v3 \
                    --sample project-GP8V3yjJXgZZ9vbKBFz74V2Q:/Bulk/Imputation/UKB imputation from genotype/ukb22828_c${chr}_b0_v3
                    --extract snpsPassInfo_0.7_${chr}.txt \
                    --remove ids_to_exclude.txt \
                    --geno 0.02 --maf 0.005 --hwe 1*10-6 \
                    --make-bed --out ${chr}_qcd" 

   dx run swiss-army-knife \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/Bulk/Imputation/UKB imputation from genotype/ukb22828_c${chr}_b0_v3.bgen" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/Bulk/Imputation/UKB imputation from genotype/ukb22828_c${chr}_b0_v3.sample" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/snpsPassInfo_0.7_${chr}.txt \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/ids_to_exclude.txt \
      -icmd"${plink_qc}" \
      --instance-type="mem3_ssd2_v2_x16" # Set node I want to use
done
```
#### Step 4: Merge QC'd plink files

First, create a file listing filesets to merge in bash, and upload this to dnanexus using dx upload.
```bash
for i in {1..22}; do echo ${i}"_qcd.bed" $i"_qcd.bim" $i"_qcd.fam"; done > allfiles.txt
```
Next, run the folling from the dx command line client:
```bash
merge="plink --merge-list allfiles.txt \
              --out qcd_merged"

dx run swiss-army-knife \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/allfiles.txt" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/1_chronotype_qcd.bed" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/1_chronotype_qcd.fam" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/1_chronotype_qcd.bim" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/2_chronotype_qcd.bed" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/2_chronotype_qcd.fam" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/2_chronotype_qcd.bim" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/3_chronotype_qcd.bed" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/3_chronotype_qcd.fam" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/3_chronotype_qcd.bim" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/4_chronotype_qcd.bed" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/4_chronotype_qcd.fam" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/4_chronotype_qcd.bim" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/5_chronotype_qcd.bed" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/5_chronotype_qcd.fam" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/5_chronotype_qcd.bim" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/6_chronotype_qcd.bed" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/6_chronotype_qcd.fam" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/6_chronotype_qcd.bim" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/7_chronotype_qcd.bed" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/7_chronotype_qcd.fam" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/7_chronotype_qcd.bim" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/8_chronotype_qcd.bed" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/8_chronotype_qcd.fam" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/8_chronotype_qcd.bim" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/9_chronotype_qcd.bed" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/9_chronotype_qcd.fam" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/9_chronotype_qcd.bim" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/10_chronotype_qcd.bed" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/10_chronotype_qcd.fam" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/10_chronotype_qcd.bim" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/11_chronotype_qcd.bed" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/11_chronotype_qcd.fam" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/11_chronotype_qcd.bim" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/12_chronotype_qcd.bed" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/12_chronotype_qcd.fam" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/12_chronotype_qcd.bim" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/13_chronotype_qcd.bed" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/13_chronotype_qcd.fam" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/13_chronotype_qcd.bim" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/14_chronotype_qcd.bed" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/14_chronotype_qcd.fam" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/14_chronotype_qcd.bim" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/15_chronotype_qcd.bed" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/15_chronotype_qcd.fam" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/15_chronotype_qcd.bim" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/16_chronotype_qcd.bed" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/16_chronotype_qcd.fam" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/16_chronotype_qcd.bim" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/17_chronotype_qcd.bed" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/17_chronotype_qcd.fam" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/17_chronotype_qcd.bim" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/18_chronotype_qcd.bed" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/18_chronotype_qcd.fam" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/18_chronotype_qcd.bim" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/19_chronotype_qcd.bed" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/19_chronotype_qcd.fam" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/19_chronotype_qcd.bim" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/20_chronotype_qcd.bed" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/20_chronotype_qcd.fam" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/20_chronotype_qcd.bim" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/21_chronotype_qcd.bed" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/21_chronotype_qcd.fam" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/21_chronotype_qcd.bim" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/22_chronotype_qcd.bed" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/22_chronotype_qcd.fam" \
      -iin="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/22_chronotype_qcd.bim" \
      -icmd="${merge}" \
      --instance-type="mem2_ssd2_v2_x32" # Set node I want to use
```
#### Step 5: Create lists of samples with schizophrenia and bipolar to be removed for the sczizophrenia and bipolar polygeinc score analyses, respectively.
```bash
### samples with scz

grep "F20" icd10_diagnoses_participant.tsv | awk '{print $1}' > samples_icd10_scz.txt

awk '$2 == 2 || $3 == 2 || $4 == 2 || $5 == 2 || $6 == 2 || $7 == 2 || $8 == 2 || $9 == 2 || $10 == 2 || $11 == 2 || $12 == 2 || $13 == 2 || $14 == 2 {print $1}' MentalHealthDiagnosedProfessional.tsv >> samples_icd10_scz.txt

sort samples_icd10_scz.txt | uniq > samples_scz.txt # These are the samples I want to remove

### samples with bp

rep -E "F31" icd10_diagnoses_participant.tsv | awk '{print $1}' > samples_icd10_bp.txt
awk '$2 == 10 || $3 == 10 || $4 == 10 || $5 == 10 || $6 == 10 || $7 == 10 || $8 == 10 || $9 == 10 || $10 == 10 || $11 == 10 || $12 == 10 || $13 == 10 || $14 == 10' MentalHealthDiagnosedProfessional.tsv | awk '{prink $1}' >> samples_icd10_bp.txt
sort samples_icd10_bp.txt | uniq > samples_bp.txt
```

## Preparation of discovery data
### For genome-wide analyses using [SBayesRC](https://www.biorxiv.org/content/10.1101/2022.10.12.510418v1)


### For pathway-based analyses using [PRSet](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010624) 

## Running PRSet
PRSet was run from the DNAnexus command line, using the app, prset_rsamples_ncovar, which was modified from the app, PRSice, already available on UKB RAP. 

Below is an example of running PRSet, using chronotype as the target data and bipolar summary statistics as the discovery data.

```bash
dx run project-GP8V3yjJXgZZ9vbKBFz74V2Q:/my_apps/prset_rsamples_ncovar --instance-type mem3_ssd1_v2_x32 -ibase_assoc="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/prset/input_files/bp_files/daner_bip_pgc3_nm_noukbiobank.filt0.5.assoc" -iplink_bed="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/chrono_merged_filt.bed" -iplink_bim="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/chrono_merged_filt.bim" -iplink_fam="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/genotype_qc/chrono_merged_filt.fam" -igtf="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/prset/input_files/Homo_sapiens.GRCh37.87.gtf" -imsigdb="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/prset/input_files/pathways502500.txt" -iextract_snps="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/prset/input_files/bp_files/snps2keep_filt0.5.bim" -iremove_samples="project-GP8V3yjJXgZZ9vbKBFz74V2Q:/prset/input_files/scz_files/samples_bp.txt" -iextra_options="--wind-3 10k --wind-5 35k --proxy 0.8 --A1 A1 --A2 A2 --pvalue P --clump-r2 0.2 --stat OR --snp SNP --base-info INFO:0.7 --set-perm 5000" --priority low
```



