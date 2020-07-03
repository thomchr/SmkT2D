library(TwoSampleMR)
library(dplyr)
library(psych)
library(MVMR)
library(MendelianRandomization)

# load exposure smoking data
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

# load exposure bmi data
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


# load outcome t2d data
pad <- read.table("/Users/dr/Desktop/MR.SMK-T2D/Mahajan.NatGenet2018b.T2D.European.hg19rsid.txt", header=T)
pad$phenotype = "T2D"
pad_data <- format_data(pad, type="outcome", snps=pad$rsid, header=T,
                        snp_col="rsid",
                        beta_col="Beta", se_col="SE",
                        eaf_col="EAF", effect_allele_col="EA",
                        other_allele_col="NEA", pval_col="Pvalue",
                        chr_col = "Chr", pos_col ="Pos")


rm(pad)
rm(bmi)
rm(smk)

# not enough memory to directly run harmonization. Reduce the data first
smk_index <- (smk_data$chr.exposure %in% bmi_data$chr.exposure) & (smk_data$pos.exposure %in% bmi_data$pos.exposure)
smk_red <- smk_data[smk_index,]
nrow(smk_red) # 2398164

bmi_index <- (bmi_data$chr.exposure %in% smk_data$chr.exposure) & (bmi_data$pos.exposure %in% smk_data$pos.exposure)
bmi_red <- bmi_data[bmi_index,]
nrow(bmi_red) # 2332179

# harmonising by T2D data
smk_harm <- harmonise_data(exposure_dat = smk_red, outcome_dat = pad_data)
nrow(smk_harm) # 2373159

bmi_harm <- harmonise_data(exposure_dat = bmi_red, outcome_dat = pad_data)
nrow(bmi_harm) # 2317364

rm(smk_red)
rm(bmi_red)
rm(smk_data)
rm(bmi_data)

smk_harm$key <- paste(smk_harm$chr.exposure, smk_harm$pos.exposure, sep = ":")
bmi_harm$key <- paste(bmi_harm$chr.exposure, bmi_harm$pos.exposure, sep = ":")

# merge, x is smoking and y is bmi
merge <- merge(smk_harm, bmi_harm, by.x = "key", by.y = "key", all.x = FALSE, all.y = FALSE)
nrow(merge) # 2317213
merge.sig <- merge[(merge$pval.exposure.x <= 5*10^-8) | (merge$pval.exposure.y <= 5*10^-8),]
nrow(merge.sig) # 44472

rm(merge)

# clumping
# get rid of the ".x" in order to run clump_data
d <- colnames(merge.sig) 
cat(paste(d, sep = '', collapse = '\",\"'))
# ****** this line may need to be changed if using different smk and/or bmi data.
colnames(merge.sig) <- c("key","SNP","effect_allele.exposure","other_allele.exposure",
                         "effect_allele.outcome","other_allele.outcome","beta.exposure",
                         "beta.outcome","eaf.exposure","eaf.outcome","remove","palindromic",
                         "ambiguous","id.outcome","chr.outcome","pos.outcome","se.outcome",
                         "pval.outcome","outcome","mr_keep.outcome","pval_origin.outcome",
                         "chr.exposure","pos.exposure","se.exposure","pval.exposure","exposure",
                         "mr_keep.exposure","pval_origin.exposure","id.exposure","action","mr_keep",
                         "samplesize.outcome",
                         "SNP.y","effect_allele.exposure.y","other_allele.exposure.y",
                         "effect_allele.outcome.y","other_allele.outcome.y","beta.exposure.y",
                         "beta.outcome.y","eaf.exposure.y","eaf.outcome.y","remove.y","palindromic.y",
                         "ambiguous.y","id.outcome.y","chr.outcome.y","pos.outcome.y","se.outcome.y",
                         "pval.outcome.y","outcome.y","mr_keep.outcome.y","pval_origin.outcome.y",
                         "chr.exposure.y","pos.exposure.y","se.exposure.y","pval.exposure.y",
                         "samplesize.exposure","exposure.y","mr_keep.exposure.y","pval_origin.exposure.y",
                         "id.exposure.y","action.y","mr_keep.y","samplesize.outcome.y")

clumped <- clump_data(merge.sig, clump_r2 = 0.01, clump_kb = 250)
nrow(clumped) # 1439

# flip the beta signs
inconsistant <- clumped$beta.outcome != clumped$beta.outcome.y
clumped[inconsistant,]$beta.exposure.y <- clumped[inconsistant,]$beta.exposure.y * -1



# perform mvmr
mvmr_in <- format_mvmr(BXGs = clumped[,c("beta.exposure", "beta.exposure.y")], 
                       BYG = clumped$beta.outcome, 
                       seBXGs = clumped[,c("se.exposure", "se.exposure.y")],
                       seBYG = clumped$se.outcome,
                       RSID = clumped$SNP)


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
                   "effect allele frequency", "lifetime smoking effect size",
                   "lifetime smoking standard error", "BMI effect size", "BMI standard error",
                   "T2D effect size", "T2D standard error")

write.csv(int, "/Users/dr/Desktop/MR.SMK-T2D/better_instrumnet_mvmr/instruments/combined_IV.LtSmk_bmi->t2d.instrument.csv", row.names = FALSE)






