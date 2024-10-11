# Load required libraries
library(dplyr)

# Read in geneset
mrna_genes <- read.table("rna_splice_minor.txt", header=F)

# Read in GTF file
gtf_genes <- read.table("gtf_genes_corord.txt", header=F) # edited this in bash

# add 35kb to start of gene
gtf_genes <- gtf_genes %>%
  mutate(
    V2 = ifelse(V4 == "+", V2 + 35000, V2),
    V3 = ifelse(V4 == "-", V3 + 35000, V3),
    V3 = ifelse(V4 == "-", V3 + 10000, V3)
  )

# Subset gtf file to only include genes in geneset
gtf_mrna <- semi_join(gtf_genes, mrna_genes, by = c("V5" = "V1"))


# add column names
colnames(gtf_mrna) <- c("CHR", "start", "stop", "strand", "gene")

# read in GWAS summary statistics
bp_ss <- read.table("~/pathway_proj/bp_sumstats/bp_ss.txt", header=T)

# Create function to subset gwas summary statistics based on whether the basepair values fall within the start and stop ranges in gtf_mrna for matching CHR values, 

is_in_range <- function(chr, bp, df1) {
  df1_subset <- df1 %>% filter(CHR == chr)
  any(df1_subset$start <= bp & bp <= df1_subset$stop)
}

# Apply this function to create bp_mrna
bp_mrna <- bp_ss %>%
  rowwise() %>%
  filter(is_in_range(CHR, BP, gtf_mrna)) %>%
  ungroup()
# convert from tibble to dataframe
bp_mrna <- as.data.frame(bp_mrna)

write.table(bp_mrna, "bp_mrna.txt", row.names=F, col.names=T, quote=F)

# read in eqtl data
brain_cerebellum <- read.table("/home/lfahey/pathway_proj/corrections/eqtl/Brain-Cerebellum.txt", header=T)

# Create SNP columns in both dataframes
brain_cerebellum$SNP <- paste(brain_cerebellum$SNP_chr, brain_cerebellum$SNP_pos, sep=":")
bp_mrna$SNP <- paste(bp_mrna$CHR, bp_mrna$BP, sep=":")

# Merge two dataframes by SNP columns
input <- merge(brain_cerebellum, bp_mrna, by="SNP", all=FALSE, suffixes=c("_eqtl","gwas"))

# Filter this to only genes in gene-set
input_mrna <- semi_join(input, mrna_genes, by = c("Mapped_gene" = "V1"))
# look at all unique genes
unique(input_mrna$Mapped_gene)

# Now, filter to just one gene
input_SF3B4 <- input_mrna[input_mrna$Mapped_gene == "SF3B4",]

# Run coloc
result <- coloc.abf(
  dataset1 = list(
    pvalues = input_SF3B4$P, # have this
    type = "cc", # will be for bp
    s = 0.1,
    N=344311,
    snp = input_SF3B4$SNP
  ),
  dataset2 = list(
    pvalues = input_SF3B4$Pvalue,
    type = "quant",
    N=10000, # NEED TO FIGURE OUT IF THIS IS DEFO PART. NO.
    snp = input_SF3B4$SNP
  ),
  MAF = input_SF3B4$FRQ_A_40463 #
)

snp_res <- result$results
head(snp_res[order(-snp_res$SNP.PP.H4),])

# Regression coefficients should be used if available. - How 2 supply these? beta=
# Write a loop for above


df_unique <- input_SF3B1[!duplicated(input_SF3B1$SNP), ]
head(df_unique)
result <- coloc.abf(
  dataset1 = list(
    pvalues = df_unique$P, # have this
    type = "cc", # will be for bp
    s = 0.1,
    N=30311,
    snp = df_unique$SNP
  ),
  dataset2 = list(
    pvalues = df_unique$Pvalue,
    type = "quant",
    N=1000,
    snp = df_unique$SNP
  ),
  MAF = df_unique$FRQ_A_40463 #
)

# Extract the more likely causal variants

