### Step 1: Load libraries and read in files
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(ggplot2)

# Read in .summary files as dataframes
asd_res <- read.table("//wsl.localhost/Ubuntu/home/lfahey/pathway_proj/all_prset_outputs/chronotype_res/asd_sumstats0.5.summary", header=T)
adhd_res <- read.table("//wsl.localhost/Ubuntu/home/lfahey/pathway_proj/all_prset_outputs/ADHD_meta_Jan2022_iPSYCH1_iPSYCH2_deCODE_PGC.filt0.5.summary", header=T)
scz_res <- read.table("//wsl.localhost/Ubuntu/home/lfahey/pathway_proj/all_prset_outputs/chronotype_res/PGC3_SCZ_wave3.primary.autosome.public.v3.filt0.5.summary", header=T)
bp_res <- read.table("//wsl.localhost/Ubuntu/home/lfahey/pathway_proj/all_prset_outputs/chronotype_res/daner_bip_pgc3_nm_noukbiobank.filt0.5.summary", header=T)

asd_res_ins <- read.table("//wsl.localhost/Ubuntu/home/lfahey/pathway_proj/all_prset_outputs/insomnia/ins_merged_filt/asd_sumstats0.5.summary", header=T)
adhd_res_ins <- read.table("//wsl.localhost/Ubuntu/home/lfahey/pathway_proj/all_prset_outputs/insomnia/ins_merged_filt/ADHD_meta_Jan2022_iPSYCH1_iPSYCH2_deCODE_PGC.filt0.5.summary", header=T)
scz_res_ins <- read.table("//wsl.localhost/Ubuntu/home/lfahey/pathway_proj/all_prset_outputs/insomnia/ins_merged_filt/PGC3_SCZ_wave3.primary.autosome.public.v3.filt0.5.summary", header=T)
bp_res_ins <- read.table("//wsl.localhost/Ubuntu/home/lfahey/pathway_proj/all_prset_outputs/insomnia/ins_merged_filt/daner_bip_pgc3_nm_noukbiobank.filt0.5.summary", header=T)

### Step 2: Prep matrix of P values to graph
asd_res_p <- asd_res[,c("Set", "Competitive.P")]
adhd_res_p <- adhd_res[,c("Set", "Competitive.P")]
scz_res_p <- scz_res[,c("Set", "Competitive.P")]
bp_res_p <- bp_res[,c("Set", "Competitive.P")]

# ins

ins_asd_res_p <- asd_res_ins[,c("Set", "Competitive.P")]
ins_adhd_res_p <- adhd_res_ins[,c("Set", "Competitive.P")]
ins_scz_res_p <- scz_res_ins[,c("Set", "Competitive.P")]
ins_bp_res_p <- bp_res_ins[,c("Set", "Competitive.P")]

# Order dataframes by competitive p-values
sig_asd_res <- head(asd_res_p[order(asd_res_p$Competitive.P, decreasing = FALSE),], 4)
sig_adhd_res <- head(adhd_res_p[order(adhd_res_p$Competitive.P, decreasing = FALSE),], 4)
sig_scz_res <- head(scz_res_p[order(scz_res_p$Competitive.P, decreasing = FALSE),], 4)
sig_bp_res <- head(bp_res_p[order(bp_res_p$Competitive.P, decreasing = FALSE),], 4)

# ins
ins_sig_asd_res <- head(ins_asd_res_p[order(ins_asd_res_p$Competitive.P, decreasing = FALSE),], 4)
ins_sig_adhd_res <- head(ins_adhd_res_p[order(ins_adhd_res_p$Competitive.P, decreasing = FALSE),], 4)
ins_sig_scz_res <- head(ins_scz_res_p[order(ins_scz_res_p$Competitive.P, decreasing = FALSE),], 4)
ins_sig_bp_res <- head(ins_bp_res_p[order(ins_bp_res_p$Competitive.P, decreasing = FALSE),], 4)


# Get list of significant sets and restrict dfs to these
sig_sets_chron_ins <- c(sig_bp_res$Set, sig_scz_res$Set, sig_asd_res$Set, sig_adhd_res$Sets, ins_sig_bp_res$Set, ins_sig_scz_res$Set, ins_sig_asd_res$Set, ins_sig_adhd_res$Set)

#chrono
asd_res_p_sig <- asd_res_p[asd_res_p$Set %in% sig_sets_chron_ins,]
adhd_res_p_sig <- adhd_res_p[adhd_res_p$Set %in% sig_sets_chron_ins,]
scz_res_p_sig <- scz_res_p[scz_res_p$Set %in% sig_sets_chron_ins,]
bp_res_p_sig <- bp_res_p[bp_res_p$Set %in% sig_sets_chron_ins,]

# ins
ins_asd_res_p_sig <- ins_asd_res_p[ins_asd_res_p$Set %in% sig_sets_chron_ins,]
ins_adhd_res_p_sig <- ins_adhd_res_p[ins_adhd_res_p$Set %in% sig_sets_chron_ins,]
ins_scz_res_p_sig <- ins_scz_res_p[ins_scz_res_p$Set %in% sig_sets_chron_ins,]
ins_bp_res_p_sig <- ins_bp_res_p[ins_bp_res_p$Set %in% sig_sets_chron_ins,]

# Change column names
# chrono
colnames(asd_res_p_sig) <- c("asd.Set","asd.P")
colnames(adhd_res_p_sig) <- c("adhd.Set","adhd.P")
colnames(scz_res_p_sig) <- c("scz.Set", "scz.P")
colnames(bp_res_p_sig) <- c("bp.Set", "bp.P")

# ins
colnames(ins_asd_res_p_sig) <- c("ins.asd.Set","ins.asd.P")
colnames(ins_adhd_res_p_sig) <- c("ins.adhd.Set","ins.adhd.P")
colnames(ins_scz_res_p_sig) <- c("ins.scz.Set", "ins.scz.P")
colnames(ins_bp_res_p_sig) <- c("ins.bp.Set", "ins.bp.P")

# Bind by columns
asd_adhd_scz_res <- cbind(adhd_res_p_sig, asd_res_p_sig, scz_res_p_sig, bp_res_p_sig, ins_asd_res_p_sig, ins_adhd_res_p_sig, ins_scz_res_p_sig, ins_bp_res_p_sig)

# Subset
asd_adhd_scz__bp_P <- asd_adhd_scz_res[,c("asd.Set", "adhd.P", "asd.P", "bp.P", "scz.P", "ins.adhd.P", "ins.asd.P", "ins.bp.P", "ins.scz.P")]

# add pathways as row names
df_res <- asd_adhd_scz__bp_P[,-1]
row.names(df_res) <- asd_adhd_scz__bp_P[,1]

# Convert to matrix
mat_res_P <- data.matrix(df_res)
rownames(mat_res_P) <- c("Antigen Processing-Cross Presentation", "Signalling by NOTCH", "NLR Signalling Pathways", "Signalling by PDGF", "Scavenging of HEME from Plasma", "Resolution of Sister Chromatid Cohesion", "DNA Damage Stress Induced Senescence", "FCERI mediated NF-kB activation", "Circadian Clock", "RHO GTPases Activate Formins", "Diseases of Metabolism", "MAPK6/MAPK4 Signaling", "PTEN Regulation", "mRNA Splicing - Minor Pathway", "Processing of Capped Intron-Containing Pre-mRNA", "RNA Polymerase I Transcription", "RAB GEFs exchange GTP for GDP on RABs", "RHOD GTPase Cycle", "RAC3 GTPase Cycle", "Signalling by NOTCH4", "Potential Therapeutics for SARS", "KEAP1-NRF2 Pathway", "Nuclear Events Mediated by NRF2", "Antigen processing: Degradation")


# Shorten row names
#rownames(mat_res_P) <- sub("REACTOME_", "", rownames(mat_res_P))
#rownames(mat_res_P) <- gsub("_", " ", rownames(mat_res_P))
#rownames(mat_res_P) <- gsub("_", " ", rownames(mat_res_P))
#rownames(mat_res_P) <- gsub("NFE2L2", "NRF2", rownames(mat_res_P))


# breaks to match colours
# Lowest threshold should be statistical significance threshold
mat_breaks <- c(0, 0.0005, 0.005, 0.05, 0.1, 0.2, 0.3, 0.4, 1)

# Test colours
br_pal <- brewer.pal(11,"RdYlBu")
scales::show_col(br_pal)
colours <- c("#A50026", "#D73027", "#FDAE61", "#FEE090", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4")

# Create heatmap and save to file
heatmap <- pheatmap(mat_res_P, border_color = "grey60", labels_col = c("AUT", "ADHD", "SCZ", "BP", "AUT", "ADHD", "SCZ", "BP"), angle_col = 0, gaps_col = c(4, 4, 4), cluster_cols = FALSE, breaks = mat_breaks, show_rownames = TRUE, legend_labels = mat_breaks, color = colorRampPalette(brewer.pal(n = 11, "RdYlBu"))(11))
changeLegend(mat_breaks, colors)

ggsave(
  "heatmap_sig_short_rownames.jpeg",
  device = "jpeg",
  plot = heatmap)


  
