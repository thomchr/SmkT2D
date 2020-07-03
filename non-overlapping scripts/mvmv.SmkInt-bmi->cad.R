library(TwoSampleMR)
library(dplyr)
library(psych)
library(MVMR)
library(MendelianRandomization)

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
bmi <- read.delim("/Users/dr/Desktop/MR.SMK-T2D/non_overlap_gwas_data/BMI_Locke_2015/GIANT_2015_BMI_EUR_withhg19pos.cp.txt", header = T)
bmi$phenotype = "bmi"
bmi_data <- format_data(bmi, type="exposure", snps=bmi$SNP, header=T,
                        snp_col="SNP", phenotype_col = "phenotype",
                        beta_col="b", se_col="se",
                        pval_col="p", samplesize_col="N",
                        eaf_col="Freq1.Hapmap",
                        effect_allele_col="A1",
                        other_allele_col="A2",
                        chr_col = "chr", pos_col = "pos")

# load outcome cad data
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
rm(cad)
rm(bmi_index)
rm(smk_index)

# harmonising by T2D data
smk_harm <- harmonise_data(exposure_dat = smk_red, outcome_dat = cad_data)
bmi_harm <- harmonise_data(exposure_dat = bmi_red, outcome_dat = cad_data)

rm(smk_red)
rm(bmi_red)

smk_harm$key <- paste(smk_harm$chr.exposure, smk_harm$pos.exposure, sep = ":")
bmi_harm$key <- paste(bmi_harm$chr.exposure, bmi_harm$pos.exposure, sep = ":")

# merge, x is smoking and y is bmi
merge <- merge(smk_harm, bmi_harm, by.x = "key", by.y = "key", all.x = FALSE, all.y = FALSE)
nrow(merge) # 2435601

merge.sig <- merge[(merge$pval.exposure.x <= 5*10^-8) | (merge$pval.exposure.y <= 5*10^-8),]
nrow(merge.sig)  # 5208


# clumping
# get rid of the ".x" in order to run clump_data
d <- colnames(merge.sig) 
cat(paste(d, sep = '', collapse = '\",\"'))
colnames(merge.sig) <- c("key","SNP","effect_allele.exposure","other_allele.exposure",
                         "effect_allele.outcome","other_allele.outcome","beta.exposure",
                         "beta.outcome","eaf.exposure","eaf.outcome","remove","palindromic",
                         "ambiguous","id.outcome","chr.outcome","pos.outcome","se.outcome",
                         "pval.outcome","outcome","mr_keep.outcome","pval_origin.outcome",
                         "chr.exposure","pos.exposure","pval.exposure","se.exposure","samplesize.exposure",
                         "exposure","mr_keep.exposure","pval_origin.exposure","id.exposure",
                         "action","mr_keep","samplesize.outcome",
                         "SNP.y","effect_allele.exposure.y","other_allele.exposure.y",
                         "effect_allele.outcome.y","other_allele.outcome.y","beta.exposure.y",
                         "beta.outcome.y","eaf.exposure.y","eaf.outcome.y","remove.y","palindromic.y",
                         "ambiguous.y","id.outcome.y","chr.outcome.y","pos.outcome.y","se.outcome.y",
                         "pval.outcome.y","outcome.y","mr_keep.outcome.y","pval_origin.outcome.y",
                         "se.exposure.y","pval.exposure.y","samplesize.exposure.y",
                         "chr.exposure.y","pos.exposure.y","exposure.y","mr_keep.exposure.y",
                         "pval_origin.exposure.y","id.exposure.y","action.y","mr_keep.y","samplesize.outcome.y")

clumped <- clump_data(merge.sig, clump_r2 = 0.01, clump_kb = 250)
nrow(clumped) # 207

# flip the beta signs
inconsistant <- clumped$beta.outcome != clumped$beta.outcome.y
clumped[inconsistant,]$beta.exposure.y <- clumped[inconsistant,]$beta.exposure.y * -1

# perform mvmr
mvmr_in <- format_mvmr(BXGs = clumped[,c("beta.exposure", "beta.exposure.y")], 
                       BYG = clumped$beta.outcome, 
                       seBXGs = clumped[,c("se.exposure", "se.exposure.y")],
                       seBYG = clumped$se.outcome,
                       RSID = clumped$SNP)

# package: MVMR
# ******note: gencov is set to default 0. 
# ******may need  covariance between the effect of the genetic variants on each exposure
mvmr_out <- mvmr(mvmr_in, gencov = 0, weights = 1)



length(which(clumped$pval.exposure <= 5*10^-8))
length(which(clumped$pval.exposure.y <= 5*10^-8))

# output the instrument
int <- clumped[, c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure", "other_allele.exposure",
                       "eaf.exposure", "beta.exposure", "se.exposure", "beta.exposure.y", "se.exposure.y",
                       "beta.outcome", "se.outcome")]

colnames(int) <- c("rsid", "chromosome", "position", "effect allele", "non-effect allele",
                   "effect allele frequency", "smoking-initiation effect size",
                   "smoking-initiation standard error", "BMI effect size", "BMI standard error",
                   "CAD effect size", "CAD standard error")

write.csv(int, "/Users/dr/Desktop/MR.SMK-T2D/non_overlap_analyses/instruments/instrument.SmkInt_bmi->cad.csv", row.names = FALSE)







