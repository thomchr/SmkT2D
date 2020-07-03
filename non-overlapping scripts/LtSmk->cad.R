library(TwoSampleMR)
library(dplyr)
library(psych)
library(stringr)

## smoking -> CAD
### Load data
# load smoking data
smk <- read.table("/Users/dr/Desktop/MR.SMK-T2D/lifetime_smoking/2019.10.02.SummStats.txt", header=T)
smk$phenotype = "smoking"
smk_data <- format_data(smk, type="exposure", snps=smk$SNP, header=T,
                        snp_col="SNP", phenotype_col = "phenotype",
                        beta_col="BETA", se_col="SE",
                        pval_col="P", 
                        eaf_col="EAF",
                        effect_allele_col="EFFECT_ALLELE",
                        other_allele_col="OTHER_ALLELE",
                        chr_col = "CHR", pos_col = "BP")
smk_data <- smk_data[smk_data$pval.exposure <  5 * 10^-8,]


# load CAD data: CAD_Nikpay_2015
cad <- read.delim("/Users/dr/Desktop/MR.SMK-T2D/non_overlap_gwas_data/CAD_Nikpay_2015/cad.add.160614.website.txt", header = T)
cad_data <- format_data(cad, type="outcome", snps=cad$markername, header=T,
                        snp_col="markername", phenotype_col = "cad",
                        beta_col="beta", se_col="se_dgc",
                        pval_col="p_dgc", 
                        eaf_col="effect_allele_freq",
                        effect_allele_col="effect_allele",
                        other_allele_col="noneffect_allele",
                        chr_col = "chr", pos_col = "bp_hg19")
cad_data$outcome <- "cad"

# harmonising
smk_harm <- harmonise_data(exposure_dat = smk_data, outcome_dat = cad_data)
nrow(smk_harm)   # 10299

# clumping
smk_clumped <- clump_data(smk_harm, clump_r2 = 0.01, clump_kb = 250)
nrow(smk_clumped)   # 194

# output instrument
int <- smk_clumped[, c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure", "other_allele.exposure",
                       "eaf.exposure", "beta.exposure", "se.exposure", "beta.outcome", "se.outcome")]

colnames(int) <- c("rsid", "chromosome", "position", "effect allele", "non-effect allele",
                   "effect allele frequency", "Lifetime-smoking effect size",
                   "Lifetime-smoking standard error", "CAD effect size", "CAD standard error")

write.csv(int, "/Users/dr/Desktop/MR.SMK-T2D/non_overlap_analyses/instruments/instrument.LtSmk->cad.csv", row.names = FALSE)


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
pdf("/Users/dr/Desktop/MR.SMK-T2D/non_overlap_analyses/scatterplots/LtSmk->cad.scatter.pdf")
p1[[1]]
dev.off()

### MR Steiger test of directionality
directionality_test(smk_clumped)
mr_steiger(p_exp = smk_clumped$pval.exposure, p_out = smk_clumped$pval.outcome, 
           n_exp = smk_clumped$samplesize.exposure, n_out = smk_clumped$samplesize.outcome, 
           r_exp = smk_clumped$beta.exposure, r_out = smk_clumped$beta.outcome, 
           r_xxo = 1, r_yyo = 1)


### r2: for getting F-statistics
smk_clumped$samplesize.exposure <- 462690
smk_clumped$r2 <- 2*(smk_clumped$beta.exposure)^2 * smk_clumped$eaf.exposure * (1-smk_clumped$eaf.exposure) /
  (2*(smk_clumped$beta.exposure)^2 * smk_clumped$eaf.exposure * (1-smk_clumped$eaf.exposure) + 
     (smk_clumped$se.exposure)^2*2*smk_clumped$samplesize.exposure*smk_clumped$eaf.exposure * (1-smk_clumped$eaf.exposure))
sum(smk_clumped$r2) # 0.009059311

smk_clumped$samplesize.exposure[1:5]



