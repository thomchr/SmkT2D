library(TwoSampleMR)
library(dplyr)
library(psych)
library(stringr)

## smoking -> T2D
### Load data
# load smoking data
smk <- read.table("/Users/dr/Desktop/MR.SMK-T2D/lifetime_smoking/2019.10.02.SummStats.txt", header=T)
smk$phenotype = "smoking"
smk_data <- format_data(smk, type="outcome", snps=smk$SNP, header=T,
                        snp_col="SNP", phenotype_col = "phenotype",
                        beta_col="BETA", se_col="SE",
                        pval_col="P", 
                        eaf_col="EAF",
                        effect_allele_col="EFFECT_ALLELE",
                        other_allele_col="OTHER_ALLELE",
                        chr_col = "CHR", pos_col = "BP")
smk_data$key <- paste(smk_data$chr.outcome, smk_data$pos.outcome, sep = ":")



# load T2D data: T2D_Scott_2017
t2d <- read.delim("/Users/dr/Desktop/MR.SMK-T2D/non_overlap_gwas_data/T2D_Scott_2017/METAANALYSIS_DIAGRAM_SE1.txt", header = T)

t2d.sig <- t2d[t2d$P.value <  5 * 10^-8,]
length(which(t2d.sig$Chr.Position %in% smk_data$key == TRUE)) #1740
t2d.sig <- t2d.sig[t2d.sig$Chr.Position %in% smk_data$key,]

t2d.sig$chr <- str_split_fixed(t2d$Chr.Position, ":", 2)[1]
t2d.sig$pos <- str_split_fixed(t2d$Chr.Position, ":", 2)[2]

# assign rsid to t2d
t2d.sig <- merge(smk_data[,c("SNP", "key")], t2d.sig, 
                 by.x = "key", by.y = "Chr.Position",
                 all.x = FALSE, all.y = FALSE)
nrow(t2d.sig)    # 1740

# format t2d data
t2d_data <- format_data(t2d.sig, type="exposure", snps=t2d.sig$SNP, header=T,
                        snp_col="SNP", phenotype_col = "t2d",
                        beta_col="Effect", se_col="StdErr",
                        pval_col="P.value", samplesize_col="TotalSampleSize",
                        #eaf_col="effect_allele_freq",
                        effect_allele_col="Allele1",
                        other_allele_col="Allele2",
                        chr_col = "chr", pos_col = "pos")


# harmonising
smk_harm <- harmonise_data(exposure_dat = t2d_data, outcome_dat = smk_data)
nrow(smk_harm)   # 1740

# clumping
smk_clumped <- clump_data(smk_harm, clump_r2 = 0.01, clump_kb = 250)
nrow(smk_clumped)   # 45

# output instrument
int <- smk_clumped[, c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure", "other_allele.exposure",
                       "eaf.outcome", "beta.exposure", "se.exposure", "beta.outcome", "se.outcome")]

colnames(int) <- c("rsid", "chromosome", "position", "effect allele", "non-effect allele",
                   "effect allele frequency", "CAD effect size",
                   "CAD standard error", "Lifetime-smoking effect size", "Lifetime-smoking standard error")

write.csv(int, "/Users/dr/Desktop/MR.SMK-T2D/non_overlap_analyses/instruments/instrument.t2d->LtSmk.csv", row.names = FALSE)



### MR
smk_bmi_res <- mr(smk_clumped, method_list = c("mr_ivw", "mr_ivw_mre", "mr_ivw_fe", 
                                               "mr_egger_regression", "mr_egger_regression_bootstrap",
                                               "mr_weighted_median"))
smk_bmi_res
smk_bmi_res[,c("b", "se", "pval")]

### Heterogeneity test
het = mr_heterogeneity(smk_clumped)
het

### Horizontal pleiotropy test
hp = mr_pleiotropy_test(smk_clumped)
hp

### Satterplot
p1 <- mr_scatter_plot(smk_bmi_res, smk_clumped)
p1
pdf("/Users/dr/Desktop/MR.SMK-T2D/non_overlap_analyses/scatterplots/t2d->LtSmk.scatter.pdf")
p1[[1]]
dev.off()

### MR Steiger test of directionality
directionality_test(smk_clumped)
mr_steiger(p_exp = smk_clumped$pval.exposure, p_out = smk_clumped$pval.outcome, 
           n_exp = smk_clumped$samplesize.exposure, n_out = smk_clumped$samplesize.outcome, 
           r_exp = smk_clumped$beta.exposure, r_out = smk_clumped$beta.outcome, 
           r_xxo = 1, r_yyo = 1)


### r2: for getting F-statistics
smk_clumped$r2 <- 2*(smk_clumped$beta.exposure)^2 * smk_clumped$eaf.outcome * (1-smk_clumped$eaf.outcome) /
  (2*(smk_clumped$beta.exposure)^2 * smk_clumped$eaf.outcome * (1-smk_clumped$eaf.outcome) + 
     (smk_clumped$se.exposure)^2*2*smk_clumped$samplesize.exposure*smk_clumped$eaf.outcome * (1-smk_clumped$eaf.outcome))
sum(smk_clumped$r2) # 0.009059311

smk_clumped$samplesize.exposure[1:5]



