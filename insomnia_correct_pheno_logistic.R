#awk 'NR==FNR{a[$1];next} $1 in a {print $0}' ../pathway_proj/ukb_checking/ins_nodup_merged.fam covar_all.samples > covar_ins.samples
#awk 'FNR==NR{a[$1]=$0;next}{print a[$1]}' covar_ins.samples ../pathway_proj/ukb_checking/ins_nodup_merged.fam > covar_ins.txt
#sed -i '1s/^/IID     AGE     SEX     CENTRE  BATCH   PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9     PC10\n/' covar_ins.txt

library(dplyr)

fam <- read.table("../pathway_proj/ukb_checking/ins_nodup_merged.fam", header=F)
covar <- read.table("covar_ins.txt", header=T)

covar$PHENO <- fam$V6

covar$IID <- as.factor(covar$IID)
covar$SEX <- as.factor(covar$SEX)
covar$CENTRE <- as.factor(covar$CENTRE)
covar$BATCH <- as.factor(covar$BATCH)
covar$PHENO <- as.factor(covar$PHENO)


covar_rescale <- mutate_if(covar, is.numeric, list(~as.numeric(scale(.))))

logistic_covar <- glm(PHENO ~ SEX + BATCH + CENTRE + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + AGE, family = binomial(link = 'logit'), data=covar, na.action = na.exclude)

summary(logistic_covar)

cor_pheno <- resid(logistic_covar)

fam$V7 <- cor_pheno


write.table(fam, "corected_fam.pheno", row.names=T, col.names=F, quote=F)

#awk '{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $8}' corected_fam.pheno > corrected_fam.txt
