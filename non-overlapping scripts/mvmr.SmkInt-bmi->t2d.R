library(TwoSampleMR)
library(dplyr)
library(psych)
library(MVMR)
library(MendelianRandomization)
library(stringr)

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

# load outcome t2d data
smk_data$key <- paste(smk_data$chr.exposure, smk_data$pos.exposure, sep = ":")
t2d <- read.delim("/Users/dr/Desktop/MR.SMK-T2D/non_overlap_gwas_data/T2D_Scott_2017/METAANALYSIS_DIAGRAM_SE1.txt", header = T)
length(which(t2d$Chr.Position %in% smk_data$key == TRUE)) #9705571
t2d <- t2d[t2d$Chr.Position %in% smk_data$key,]

t2d$chr <- str_split_fixed(t2d$Chr.Position, ":", 2)[,1]
t2d$pos <- str_split_fixed(t2d$Chr.Position, ":", 2)[,2]

# assign rsid to t2d
t2d <- merge(smk_data[,c("SNP", "key")], t2d, 
             by.x = "key", by.y = "Chr.Position",
             all.x = FALSE, all.y = FALSE)
nrow(t2d)    # 9705571

# format t2d data
t2d_data <- format_data(t2d, type="outcome", snps=t2d$SNP, header=T,
                        snp_col="SNP", phenotype_col = "t2d",
                        beta_col="Effect", se_col="StdErr",
                        pval_col="P.value", samplesize_col="TotalSampleSize",
                        #eaf_col="effect_allele_freq",
                        effect_allele_col="Allele1",
                        other_allele_col="Allele2",
                        chr_col = "chr", pos_col = "pos")
t2d_data$outcome <- "t2d"


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
rm(t2d)
rm(bmi_index)
rm(smk_index)

# not enough memory to directly run harmonization. Reduce the t2d data first
t2d_index <- (t2d_data$chr.outcome %in% smk_red$chr.exposure) & (t2d_data$pos.outcome %in% smk_red$pos.exposure)
length(which(t2d_index))
t2d_red <- t2d_data[t2d_index,]

rm(t2d_data)
rm(t2d_index)


# need to further reduce
smk_index <- (smk_red$chr.exposure %in% t2d_red$chr.outcome) & (smk_red$pos.exposure %in% t2d_red$pos.outcome)
length(which(smk_index))
smk_red <- smk_red[smk_index,]
bmi_index <- (bmi_red$chr.exposure %in% t2d_red$chr.outcome) & (bmi_red$pos.exposure %in% t2d_red$pos.outcome)
length(which(bmi_index))
bmi_red <- bmi_red[bmi_index,]

# merge, x is smoking and y is bmi
smk_red$key <- paste(smk_red$chr.exposure, smk_red$pos.exposure, sep = ":")
bmi_red$key <- paste(bmi_red$chr.exposure, bmi_red$pos.exposure, sep = ":")
temp_merge <- merge(smk_red, bmi_red, by.x = "key", by.y = "key", all.x = FALSE, all.y = FALSE)
nrow(temp_merge)

temp_merge <- temp_merge[(temp_merge$pval.exposure.x <= 5*10^-8) | (temp_merge$pval.exposure.y <= 5*10^-8),]
nrow(temp_merge) # 5215

smk_sig <- temp_merge[,c("key","chr.exposure.x","pos.exposure.x","SNP.x", 
                         "other_allele.exposure.x","effect_allele.exposure.x","eaf.exposure.x","pval.exposure.x",
                         "beta.exposure.x","se.exposure.x","samplesize.exposure.x","exposure.x",
                         "mr_keep.exposure.x","pval_origin.exposure.x","id.exposure.x")]
colnames(smk_sig) <- c("key","chr.exposure","pos.exposure","SNP", 
                         "other_allele.exposure","effect_allele.exposure","eaf.exposure","pval.exposure",
                         "beta.exposure","se.exposure","samplesize.exposure","exposure",
                         "mr_keep.exposure","pval_origin.exposure","id.exposure")

bmi_sig <- temp_merge[,c("key","chr.exposure.y","pos.exposure.y","SNP.y", 
                         "other_allele.exposure.y","effect_allele.exposure.y","eaf.exposure.y","pval.exposure.y",
                         "beta.exposure.y","se.exposure.y","samplesize.exposure.y","exposure.y",
                         "mr_keep.exposure.y","pval_origin.exposure.y","id.exposure.y")]
colnames(bmi_sig) <- c("key","chr.exposure","pos.exposure","SNP", 
                       "other_allele.exposure","effect_allele.exposure","eaf.exposure","pval.exposure",
                       "beta.exposure","se.exposure","samplesize.exposure","exposure",
                       "mr_keep.exposure","pval_origin.exposure","id.exposure")

rm(bmi_red)
rm(smk_red)
rm(bmi_index)
rm(smk_index)

# harmonising by T2D data
smk_harm <- harmonise_data(exposure_dat = smk_sig, outcome_dat = t2d_red)
bmi_harm <- harmonise_data(exposure_dat = bmi_sig, outcome_dat = t2d_red)


# merge, x is smoking and y is bmi
merge <- merge(smk_harm, bmi_harm, by.x = "key", by.y = "key", all.x = FALSE, all.y = FALSE)
nrow(merge) # 5215

merge.sig <- merge


# clumping
# get rid of the ".x" in order to run clump_data
d <- colnames(merge.sig) 
cat(paste(d, sep = '', collapse = '\",\"'))
colnames(merge.sig) <- c("key","SNP","effect_allele.exposure","other_allele.exposure",
                         "effect_allele.outcome","other_allele.outcome","beta.exposure",
                         "beta.outcome","eaf.exposure","eaf.outcome","remove","palindromic",
                         "ambiguous","id.outcome","se.outcome","pval.outcome","samplesize.outcome",
                         "chr.outcome","pos.outcome","outcome","mr_keep.outcome",
                         "pval_origin.outcome","chr.exposure","pos.exposure","pval.exposure",
                         "se.exposure","samplesize.exposure","exposure","mr_keep.exposure",
                         "pval_origin.exposure","id.exposure","action","mr_keep",
                         "SNP.y","effect_allele.exposure.y","other_allele.exposure.y",
                         "effect_allele.outcome.y","other_allele.outcome.y","beta.exposure.y",
                         "beta.outcome.y","eaf.exposure.y","eaf.outcome.y","remove.y","palindromic.y",
                         "ambiguous.y","id.outcome.y","se.outcome.y","pval.outcome.y",
                         "samplesize.outcome.y","chr.outcome.y","pos.outcome.y","outcome.y",
                         "mr_keep.outcome.y","pval_origin.outcome.y","chr.exposure.y",
                         "pos.exposure.y","pval.exposure.y","se.exposure.y","samplesize.exposure.y",
                         "exposure.y","mr_keep.exposure.y","pval_origin.exposure.y","id.exposure.y",
                         "action.y","mr_keep.y")

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
                   "T2D effect size", "T2D standard error")

write.csv(int, "/Users/dr/Desktop/MR.SMK-T2D/non_overlap_analyses/instruments/instrument.SmkInt_bmi->t2d.csv", row.names = FALSE)




