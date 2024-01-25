# 1. Read in libraries
library(dplyr)
library(confintr)

# 2. Read in pgs score file (output of PRSet/PRSice/SBayesRC) and covariate file

bp_scores <- read.table("bp_chronotypePRS.profile", header=T)
covar <- read.table("/home/lfahey/ukb/covar_all.samples", header=T, fill=T)

# Covert to factors
bp_scores$IID <- as.factor(bp_scores$IID)
covar$IID <- as.factor(covar$IID)
covar$SEX <- as.factor(covar$SEX)
covar$CENTRE <- as.factor(covar$CENTRE)
covar$BATCH <- as.factor(covar$BATCH)
# include pheno for insomnia

# 3. Merge the two dataframes
bp_scores_covar <- merge(bp_scores, covar, by = "IID", all.x = FALSE)
# subset
bp_scores_covar <- bp_scores_covar[,c("IID", "PHENO", "SCORE", "AGE", "SEX", "CENTRE", "BATCH", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]

# covar_rescale <- mutate_if(covar, is.numeric, list(~as.numeric(scale(.))))

# Run a linear regression with chronotype as the dependent variable and the confounders as the independent variables
linear_covar <- lm(PHENO ~ SEX + BATCH + CENTRE + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + AGE, data=bp_scores_covar)

summary(logistic_covar)

# assign residuals as the corrected phenoytpe
bp_scores_covar$cor_pheno <- resid(linear_covar)

# Run a linear regression with the corrected phenotype as the dependent variable and the PGS as the independent variable
bp_pgs_cor <- lm(cor_pheno ~ SCORE, data = bp_scores_covar)

# Get P valaue and R2
summary(bp_pgs_cor)

# Get CIs of R2
ci_rsquared(bp_pgs_cor)
