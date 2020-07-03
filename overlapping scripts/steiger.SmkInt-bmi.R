library(TwoSampleMR)
library(dplyr)
library(psych)
library(MVMR)


# load exposure smoking data
smk <- read.table("/Users/dr/Desktop/MR.SMK-T2D/smoking_initiation_complete_snp/SmokingInitiation.txt", header=T)
smk$phenotype = "smoking"
smk_data <- format_data(smk, type="exposure", snps=smk$RSID, header=T,
                        snp_col="RSID", phenotype_col = "phenotype",
                        beta_col="BETA", se_col="SE",
                        pval_col="PVALUE", samplesize_col ="N",
                        eaf_col="AF",
                        effect_allele_col="ALT",
                        other_allele_col="REF",
                        chr_col = "CHROM", pos_col = "POS")

# load exposure bmi data
# load outcome bmi data
bmi <- read.table("/Users/dr/Desktop/MR.SMK-T2D/BMI-Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt", header=T)
bmi$phenotype = "BMI"
bmi_data <- format_data(bmi, type="exposure", snps=bmi$SNP, header=T,
                        snp_col="SNP", phenotype_col = "phenotype",
                        beta_col="BETA", se_col="SE",
                        pval_col="P", samplesize_col ="N",
                        eaf_col="Freq_Tested_Allele_in_HRS",
                        effect_allele_col="Tested_Allele",
                        other_allele_col="Other_Allele",
                        chr_col = "CHR", pos_col = "POS")


# not enough memory to directly run harmonization. Reduce the data first
smk_index <- (smk_data$chr.exposure %in% bmi_data$chr.exposure) & (smk_data$pos.exposure %in% bmi_data$pos.exposure)
smk_red <- smk_data[smk_index,]
bmi_index <- (bmi_data$chr.exposure %in% smk_data$chr.exposure) & (bmi_data$pos.exposure %in% smk_data$pos.exposure)
bmi_red <- bmi_data[bmi_index,]

# need to clear all unneeded objects, otherwise rstudio will crash
rm(bmi)
rm(bmi_data)
rm(smk)
rm(smk_data)
rm(bmi_index)
rm(smk_index)


# merge smk and bmi


head(bmi_red)
head(smk_red)

smk_red$key <- paste(smk_red$chr.exposure, smk_red$pos.exposure, sep = ":")
bmi_red$key <- paste(bmi_red$chr.exposure, bmi_red$pos.exposure, sep = ":")


merge <- merge(smk_red, bmi_red, by.x = "key", by.y = "key", all.x = FALSE, all.y = FALSE)
nrow(merge)  # 2308073


merge.sig <- merge[(merge$pval.exposure.x <= 5*10^-8) & (merge$pval.exposure.y <= 5*10^-8),]
nrow(merge.sig)  # 42445

# clump
# get rid of the ".x" in order to run clump_data
d <- colnames(merge.sig) 
cat(paste(d, sep = '', collapse = '\",\"'))
colnames(merge.sig) <- c("key","chr.exposure","pos.exposure","SNP","other_allele.exposure",
                         "effect_allele.exposure","eaf.exposure","pval.exposure",
                         "beta.exposure","se.exposure","samplesize.exposure","exposure",
                         "mr_keep.exposure","pval_origin.exposure","id.exposure",
                         "chr.exposure.y","pos.exposure.y","SNP.y","effect_allele.exposure.y",
                         "other_allele.exposure.y","eaf.exposure.y","beta.exposure.y",
                         "se.exposure.y","pval.exposure.y","samplesize.exposure.y",
                         "exposure.y","mr_keep.exposure.y","pval_origin.exposure.y","id.exposure.y")

clumped <- clump_data(merge.sig, clump_r2 = 0.01, clump_kb = 250)
nrow(clumped) # 1411


### MR Steiger test of directionality

# smk->bmi
mr_steiger(p_exp = clumped$pval.exposure, p_out = clumped$pval.exposure.y, 
           n_exp = clumped$samplesize.exposure, n_out = clumped$samplesize.exposure.y, 
           r_exp = clumped$beta.exposure, r_out = clumped$beta.exposure.y, 
           r_xxo = 0.8, r_yyo = 0.8)

# bmi->smk
mr_steiger(p_exp = clumped$pval.exposure.y, p_out = clumped$pval.exposure, 
           n_exp = clumped$samplesize.exposure.y, n_out = clumped$samplesize.exposure, 
           r_exp = clumped$beta.exposure.y, r_out = clumped$beta.exposure, 
           r_xxo = 0.8, r_yyo = 0.8)




