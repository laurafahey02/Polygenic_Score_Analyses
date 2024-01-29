library(ggplot2)

# Read in file
global_r2 <- read.table("//wsl.localhost/Ubuntu/home/lfahey/pathway_proj/results/posthoc/global_r2.txt", header=T)

### Chronotype

chronotype_r2 <- global_r2[global_r2$BASE_PHENOTYPE == "Chronotype",]

# Add levels
chronotype_r2$GWAS_PHENOTYPE <- factor(chronotype_r2$GWAS_PHENOTYPE, levels = chronotype_r2$GWAS_PHENOTYPE[order(chronotype_r2$R2)])


chrono <- ggplot(chronotype_r2, aes(x = GWAS_PHENOTYPE, y=R2, fill = BASE_PHENOTYPE)) +
  geom_bar(stat="identity") +
  scale_fill_manual(name = "Base Phenotype", values = rev(c("#76B7B2", "#B07AA1")), guide = guide_legend()) +
  geom_errorbar(aes(ymin=LOWER, ymax=UPPER,x=GWAS_PHENOTYPE),width=0.4,cex=0.5) +
  xlab('GWAS Phenotype')+ ylab(expression("R"^2*~"(95% Confidence Interval)")) +
  ggtitle("Chronotype") +
  # facet_wrap makes panels # scales = "free_y", means x-axis scale is fixed, but y axis can vary 
  coord_flip() + # flips plot axes_y
  theme_bw() +
  ylim(0.000, 0.00184) +
  theme(plot.title=element_text(size=10,face="bold",hjust = 0.5),
        axis.text.y=element_text(colour = "black", size=9), # no ticks on y-axis
        #axis.ticks.y=element_blank(), # no labels on y-axis
        axis.text.x=element_text(colour = "black", size=9), # face="bold" # change format of x-axis values
        axis.title=element_text(colour = "black", size=10), # formatting of x-axis label
        #strip.text.x = element_text(face="bold", size=12),
        strip.background = element_rect(colour="black", fill="white"))

### Insomnia

insomnia_r2 <- global_r2[global_r2$BASE_PHENOTYPE == "Insomnia",]

insomnia_r2$GWAS_PHENOTYPE <- factor(insomnia_r2$GWAS_PHENOTYPE, levels = insomnia_r2$GWAS_PHENOTYPE[order(insomnia_r2$R2)])
ins <- ggplot(insomnia_r2, aes(x = GWAS_PHENOTYPE, y=R2, fill = BASE_PHENOTYPE)) +
  geom_bar(stat="identity") +
  scale_fill_manual(name = "Base Phenotype", values = rev(c("#B07AA1", "#76B7B2")), guide = guide_legend()) +
  geom_errorbar(aes(ymin=LOWER, ymax=UPPER,x=GWAS_PHENOTYPE),width=0.1,cex=0.5) +
  xlab('GWAS Phenotype') + ylab(expression("R"^2*~"(95% Confidence Interval)")) +
  ggtitle("Insomnia") +
  coord_flip() + # flips plot axes_y
  theme_bw() +
  ylim(0.000, 0.00184) +
  theme(plot.title=element_text(size=10,face="bold",hjust = 0.5),
        axis.text.y=element_text(colour = "black", size=9), # no ticks on y-axis
        #axis.ticks.y=element_blank(), # no labels on y-axis
        axis.text.x=element_text(colour = "black", size=9), # face="bold" # change format of x-axis values
        axis.title=element_text(colour = "black", size=10), # formatting of x-axis label
        #strip.text.x = element_text(face="bold", size=12),
        strip.background = element_rect(colour="black", fill="white"))


# Save plots

ggsave("chronotype_plot.jpeg", chrono)    
