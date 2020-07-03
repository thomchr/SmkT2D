library(TwoSampleMR)
library(dplyr)
library(psych)
library(stringr)

#### get the intersection of the four datasets ####
# load lifetime smoking
smk_int <- read.table("/Users/dr/Desktop/MR.SMK-T2D/lifetime_smoking/2019.10.02.SummStats.txt", header=T)
smk_int$SNPID <- paste(smk_int$CHR, smk_int$BP, sep = ":")
nrow(smk_int) # 7683352
# smk$phenotype = "smoking"
# smk_data <- format_data(smk, type="exposure", snps=smk$SNP, header=T,
#                         snp_col="SNP", phenotype_col = "phenotype",
#                         beta_col="BETA", se_col="SE",
#                         pval_col="P", 
#                         eaf_col="EAF",
#                         effect_allele_col="EFFECT_ALLELE",
#                         other_allele_col="OTHER_ALLELE",
#                         chr_col = "CHR", pos_col = "BP")


# merge with smoking initiation
# load smoking initiation
lt_smk <- read.table("/Users/dr/Desktop/MR.SMK-T2D/smoking_initiation_complete_snp/SmokingInitiation.txt", header=T)
lt_smk$SNPID <- paste(lt_smk$CHROM, lt_smk$POS, sep = ":")
nrow(lt_smk) # 11802365
# smk$phenotype = "smoking"
# smk_data <- format_data(smk, type="outcome", snps=smk$RSID, header=T,
#                         snp_col="RSID", phenotype_col = "phenotype",
#                         beta_col="BETA", se_col="SE",
#                         pval_col="PVALUE", samplesize_col ="N",
#                         eaf_col="AF",
#                         effect_allele_col="ALT",
#                         other_allele_col="REF",
#                         chr_col = "CHROM", pos_col = "POS")

merge1 <- merge(smk_int[c("SNPID", "SNP", "EFFECT_ALLELE", "OTHER_ALLELE", "EAF","BETA", "SE", "P")], 
                lt_smk[c("SNPID", "ALT", "REF", "AF", "BETA", "SE", "PVALUE")], by.x = "SNPID", by.y = "SNPID", 
                all.x = FALSE, all.y = FALSE)
nrow(merge1) # 7539558
colnames(merge1) <- c("SNPID", "RSID", "effect.a", "other.a", "eaf.a", "beta.a", "se.a", "pval.a",
                      "effect.b", "other.b", "eaf.b", "beta.b", "se.b", "pval.b")
rm(smk_int)
rm(lt_smk)

# merge with BMI1
# load BMI: Locke_2015
bmi_locke <- read.delim("/Users/dr/Desktop/MR.SMK-T2D/non_overlap_gwas_data/BMI_Locke_2015/GIANT_2015_BMI_EUR_withhg19pos.cp.txt", header = T)
bmi_locke$SNPID <- paste(bmi_locke$chr, bmi_locke$pos, sep = ":")
nrow(bmi_locke) # 2553130
# bmi_data <- format_data(bmi, type="exposure", snps=bmi$SNP, header=T,
#                         snp_col="SNP", phenotype_col = "bmi",
#                         beta_col="b", se_col="se",
#                         pval_col="p", samplesize_col="N",
#                         eaf_col="Freq1.Hapmap",
#                         effect_allele_col="A1",
#                         other_allele_col="A2",
#                         chr_col = "chr", pos_col = "pos")
# bmi_data <- bmi_data[bmi_data$pval.exposure <  5 * 10^-8,]
merge2 <- merge(merge1, bmi_locke[c("SNPID", "A1", "A2", "Freq1.Hapmap", "b", "se", "p")], 
                by.x = "SNPID", by.y = "SNPID", 
                all.x = FALSE, all.y = FALSE)
nrow(merge2) # 2429318
colnames(merge2) <- c("SNPID", "RSID", "effect.a", "other.a", "eaf.a", "beta.a", "se.a", "pval.a",
                      "effect.b", "other.b", "eaf.b", "beta.b", "se.b", "pval.b",
                      "effect.c", "other.c", "eaf.c", "beta.c", "se.c", "pval.c")
rm(bmi_locke)


# merge with BMI2
# load BMI: Yengo
bmi_yengo <- read.table("/Users/dr/Desktop/MR.SMK-T2D/BMI-Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt", header=T)
bmi_yengo$SNPID <- paste(bmi_yengo$CHR, bmi_yengo$POS, sep = ":")
nrow(bmi_yengo) # 2336269
# bmi$phenotype = "BMI"
# bmi <- format_data(bmi, type="outcome", snps=bmi$SNP, header=T, 
#                    snp_col="SNP", phenotype_col = "phenotype",
#                    beta_col="BETA", se_col="SE",
#                    pval_col="P", samplesize_col ="N",
#                    eaf_col="Freq_Tested_Allele_in_HRS",
#                    effect_allele_col="Tested_Allele", 
#                    other_allele_col="Other_Allele")
merge3 <- merge(merge2, bmi_yengo[c("SNPID", "Tested_Allele", "Other_Allele", "Freq_Tested_Allele_in_HRS", "BETA", "SE", "P")], 
                by.x = "SNPID", by.y = "SNPID", 
                all.x = FALSE, all.y = FALSE)
nrow(merge3) # 2303817
colnames(merge3) <- c("SNPID", "RSID", "effect.a", "other.a", "eaf.a", "beta.a", "se.a", "pval.a",
                      "effect.b", "other.b", "eaf.b", "beta.b", "se.b", "pval.b",
                      "effect.c", "other.c", "eaf.c", "beta.c", "se.c", "pval.c",
                      "effect.d", "other.d", "eaf.d", "beta.d", "se.d", "pval.d")

rm(bmi_yengo)

merge3$chr <- str_split_fixed(merge3$SNPID, ":", 2)[,1]
merge3$pos <- str_split_fixed(merge3$SNPID, ":", 2)[,2]

# harmonise all the data set to smoking initiation
smk_int_data <- format_data(merge3, type="outcome", snps=merge3$RSID, header=T,
                        snp_col="RSID", beta_col="beta.a", se_col="se.a", pval_col="pval.a",
                        eaf_col="eaf.a", effect_allele_col="effect.a", other_allele_col="other.a",
                        chr_col = "chr", pos_col = "pos")

# lifetime smoking
lt_smk_data <- format_data(merge3, type="exposure", snps=merge3$RSID, header=T,
                           snp_col="RSID", beta_col="beta.b", se_col="se.b", pval_col="pval.b",
                           eaf_col="eaf.b", effect_allele_col="effect.b", other_allele_col="other.b",
                           chr_col = "chr", pos_col = "pos")


# BMI: Locke_2015
bmi_locke_data <- format_data(merge3, type="exposure", snps=merge3$RSID, header=T,
                           snp_col="RSID", beta_col="beta.c", se_col="se.c", pval_col="pval.c",
                           eaf_col="eaf.c", effect_allele_col="effect.c", other_allele_col="other.c",
                           chr_col = "chr", pos_col = "pos")

# BMI: Yengo
bmi_yengo_data <- format_data(merge3, type="exposure", snps=merge3$RSID, header=T,
                           snp_col="RSID", beta_col="beta.d", se_col="se.d", pval_col="pval.d",
                           eaf_col="eaf.d", effect_allele_col="effect.d", other_allele_col="other.d",
                           chr_col = "chr", pos_col = "pos")


harm1 <- harmonise_data(exposure_dat = lt_smk_data, outcome_dat = smk_int_data)
harm2 <- harmonise_data(exposure_dat = bmi_locke_data, outcome_dat = smk_int_data)
harm3 <- harmonise_data(exposure_dat = bmi_yengo_data, outcome_dat = smk_int_data)

colnames(harm2) <- paste(colnames(harm2), "x", sep = ".")
colnames(harm3) <- paste(colnames(harm3), "y", sep = ".")

merge4 <- merge(harm1, harm2, 
                by.x = "SNP", by.y = "SNP.x", 
                all.x = FALSE, all.y = FALSE)
nrow(merge4) # 2303790

merge5 <- merge(merge4, harm3, 
                by.x = "SNP", by.y = "SNP.y", 
                all.x = FALSE, all.y = FALSE)
nrow(merge5) # 2303790
merge5$SNPID <- paste(merge5$chr.outcome, merge5$pos.outcome, sep = ":")

# flip the beta signs
inconsistant1 <- merge5$beta.outcome != merge5$beta.outcome.x
merge5[inconsistant1,]$beta.exposure.x <- merge5[inconsistant1,]$beta.exposure.x  * -1

inconsistant2 <- merge5$beta.outcome != merge5$beta.outcome.y
merge5[inconsistant2,]$beta.exposure.y <- merge5[inconsistant2,]$beta.exposure.y  * -1

nrow(merge5) # 2303790

#### build instrument for smoking initiation ####
# keep only SNPs that are significantly associated with smoking initiation

#smk_int.sig <- merge5[merge5$pval.outcome <= 5*10^-8,]
smk_int.sig <- merge5[(merge5$pval.outcome <= 5*10^-8) | (merge5$pval.exposure <= 5*10^-8) |
                        (merge5$pval.exposure.x <= 5*10^-8) | (merge5$pval.exposure.y <= 5*10^-8),]
nrow(smk_int.sig) # 44565

# clumping
smk_int.clumped <- clump_data(smk_int.sig, clump_r2 = 0.01, clump_kb = 250)
nrow(smk_int.clumped) # 1480

# take important columns

#### clumped instrument
smk_int.inst <- smk_int.clumped[,c("SNPID", "chr.outcome", "pos.outcome", "SNP",
                                   "effect_allele.outcome", "other_allele.exposure", "eaf.outcome", 
                                   "beta.exposure.y", "se.exposure.y", "pval.exposure.y",
                                   "beta.exposure.x", "se.exposure.x", "pval.exposure.x",
                                   "beta.exposure", "se.exposure", "pval.exposure",
                                   "beta.outcome", "se.outcome", "pval.outcome")]

colnames(smk_int.inst) <- c("SNPID", "CHR", "POS", "RSID", "EFFECT_ALLELE", "OTHER_ALLELE", "EFFECT_ALLELE_FREQ",
                            "BETA_BMI_YENGO", "SE_BMI_YENGO", "P_BMI_YENGO",
                            "BETA_BMI_LOCKE", "SE_BMI_LOCKE", "P_BMI_LOCKE",
                            "BETA_SMOKING_LIFETIME", "SE_SMOKING_LIFETIME", "P_SMOKING_LIFETIME",
                            "BETA_SMOKING_INIT", "SE_SMOKING_INIT", "P_SMOKING_INIT")
write.csv(smk_int.inst, "/Users/dr/Desktop/MR.SMK-T2D/individual_level_instrument/instrument.clumped.csv", row.names = FALSE)


#### unclumped instrument
unclumped <- smk_int.sig[,c("SNPID", "chr.outcome", "pos.outcome", "SNP",
                                   "effect_allele.outcome", "other_allele.exposure", "eaf.outcome", 
                                   "beta.exposure.y", "se.exposure.y", "pval.exposure.y",
                                   "beta.exposure.x", "se.exposure.x", "pval.exposure.x",
                                   "beta.exposure", "se.exposure", "pval.exposure",
                                   "beta.outcome", "se.outcome", "pval.outcome")]

colnames(unclumped) <- c("SNPID", "CHR", "POS", "RSID", "EFFECT_ALLELE", "OTHER_ALLELE", "EFFECT_ALLELE_FREQ",
                            "BETA_BMI_YENGO", "SE_BMI_YENGO", "P_BMI_YENGO",
                            "BETA_BMI_LOCKE", "SE_BMI_LOCKE", "P_BMI_LOCKE",
                            "BETA_SMOKING_LIFETIME", "SE_SMOKING_LIFETIME", "P_SMOKING_LIFETIME",
                            "BETA_SMOKING_INIT", "SE_SMOKING_INIT", "P_SMOKING_INIT")
write.csv(unclumped, "/Users/dr/Desktop/MR.SMK-T2D/individual_level_instrument/instrument.unclumped.csv", row.names = FALSE)

nrow(unclumped)


