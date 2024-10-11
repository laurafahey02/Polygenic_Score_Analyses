# 1. Read in libraries
library(dplyr)
library(confintr)

# 2. Read in files
fam <- read.table("//wsl.localhost/Ubuntu/home/lfahey/pathway_proj/corrections/safe_fam/chrono_nodup_merged.fam", header=F)
covar <- read.table("//wsl.localhost/Ubuntu/home/lfahey/ukb/covar.txt", header=T)

fam_cols <- fam[,c("V1", "V6")]
colnames(fam_cols) <- c("IID", "PHENO")
fam_covar <- merge(fam_cols, covar, by = "IID")

fam_covar$IID <- as.factor(fam_covar$IID)
fam_covar$SEX <- as.factor(fam_covar$SEX)
fam_covar$CENTRE <- as.factor(fam_covar$CENTRE)
fam_covar$BATCH <- as.factor(fam_covar$BATCH)
fam_covar$PHENO <- as.factor(fam_covar$PHENO)

# Rescale numeric columns
fam_covar_rescale <- mutate_if(fam_covar, is.numeric, list(~as.numeric(scale(.))))
model <- clm(PHENO ~ SEX + BATCH + CENTRE + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + AGE, data=fam_covar_rescale, link = "logit")

predicted_probs <- predict(model, type = "prob")  # Linear predictors

# Calculate residuals manually: 
# Residuals can be computed as the difference between the observed response categories and the predicted categories.
observed <- as.numeric(fam_covar_rescale$PHENO)
predicted <- predicted_probs
residuals_clm <- observed - predicted$fit

# Write new fam file

fam$PHENO_COR <- residuals_clm

write.table(fam, "//wsl.localhost/Ubuntu/home/lfahey/pathway_proj/corrections/chrono_nodup_merged.fam", col.names=F, row.names=F, quote=F)
