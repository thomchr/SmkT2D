library(MVMR)
library(TwoSampleMR)
library(dplyr)

setwd("/Users/dr/Desktop/MR-SMK.T2D")

#### smoking_initiation ####
# load initiation data
smk_int <-read_exposure_data("smoking_initiation/smoking_initiation.csv", sep=",")



# already loaded
# load T2D data
#pad <- read.table("/Users/dr/Desktop/MR-SMK.T2D/Mahajan.NatGenet2018b.T2D.European.hg19rsid.txt", header=T)
#pad_data <- format_data(pad, type="outcome", snps=pad$rsid, header=T, 
#                        snp_col="rsid",
#                        beta_col="Beta", se_col="SE",
#                        eaf_col="EAF", effect_allele_col="EA",
#                        other_allele_col="NEA", pval_col="Pvalue")

# bmi data
bmi <- read.table("/Users/dr/Desktop/MR-SMK.T2D/BMI-Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt", header=T)
bmi$phenotype = "BMI"
bmi <- format_data(bmi, type="exposure", snps=bmi$SNP, header=T, 
                        snp_col="SNP", phenotype_col = "phenotype",
                        beta_col="BETA", se_col="SE",
                        pval_col="P", samplesize_col ="N",
                        eaf_col="Freq_Tested_Allele_in_HRS",
                        effect_allele_col="Tested_Allele", 
                        other_allele_col="Other_Allele")
bmi_sub <- bmi[bmi$SNP %in% smk_int$SNP,]


# LD clumping
smk_int_clumped <- clump_data(smk_int, clump_r2 = 0.01, clump_kb = 250)
bmi_clumped <- clump_data(bmi_sub, clump_r2 = 0.01, clump_kb = 250)

# harmonising
smk_int_clumped <- harmonise_data(exposure_dat = smk_int_clumped, outcome_dat = pad_data)
bmi_clumped <- harmonise_data(exposure_dat = bmi_clumped, outcome_dat = pad_data)


# change colnames for merging
colnames(smk_int_clumped)[c(6, 20)] <- c("beta.smk", "se.smk")
colnames(bmi_clumped)[c(6, 19)] <- c("beta.bmi", "se.bmi")


# merge smk, bmi, and t2d
temp <- merge(select(smk_int_clumped, SNP, beta.smk, se.smk), 
              select(pad_data, SNP, beta.outcome, se.outcome), by="SNP")
nrow(temp)
combined <- merge(select(temp, SNP, beta.smk, se.smk, beta.outcome, se.outcome), 
              select(bmi_clumped, SNP, beta.bmi, se.bmi), by="SNP")
nrow(combined)

# perform mvmr
mvmr_in <- format_mvmr(BXGs = select(combined, beta.smk, beta.bmi), 
            BYG = combined$beta.outcome, 
            seBXGs = select(combined, se.smk, se.bmi),
            seBYG = combined$se.outcome,
            RSID = combined$SNP)

# ******note: gencov is set to default 0. 
# ******may need  covariance between the effect of the genetic variants on each exposure
mvmr_out <- mvmr(mvmr_in, gencov = cov(combined$beta.smk, combined$beta.bmi), weights = 1)

#### smoking cessation ####
# load initiation data
smk_ces <-read_exposure_data("smoking_cessation/smk_ces.csv", sep=",")



# already loaded
# load T2D data
#pad <- read.table("/Users/dr/Desktop/MR-SMK.T2D/Mahajan.NatGenet2018b.T2D.European.hg19rsid.txt", header=T)
#pad_data <- format_data(pad, type="outcome", snps=pad$rsid, header=T, 
#                        snp_col="rsid",
#                        beta_col="Beta", se_col="SE",
#                        eaf_col="EAF", effect_allele_col="EA",
#                        other_allele_col="NEA", pval_col="Pvalue")

# already loaded
# bmi data
bmi <- read.table("/Users/dr/Desktop/MR-SMK.T2D/BMI-Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt", header=T)
bmi$phenotype = "BMI"
bmi <- format_data(bmi, type="exposure", snps=bmi$SNP, header=T, 
                   snp_col="SNP", phenotype_col = "phenotype",
                   beta_col="BETA", se_col="SE",
                   pval_col="P", samplesize_col ="N",
                   eaf_col="Freq_Tested_Allele_in_HRS",
                   effect_allele_col="Tested_Allele", 
                   other_allele_col="Other_Allele")
bmi_sub_ces <- bmi[bmi$SNP %in% smk_ces$SNP,]


# LD clumping
smk_ces_clumped <- clump_data(smk_ces, clump_r2 = 0.01, clump_kb = 250)
bmi_ces_clumped <- clump_data(bmi_sub_ces, clump_r2 = 0.01, clump_kb = 250)

# harmonising
smk_ces_clumped <- harmonise_data(exposure_dat = smk_ces_clumped, outcome_dat = pad_data)
bmi_ces_clumped <- harmonise_data(exposure_dat = bmi_ces_clumped, outcome_dat = pad_data)


# change colnames for merging
colnames(smk_ces_clumped)[c(6, 20)] <- c("beta.smk", "se.smk")
colnames(bmi_ces_clumped)[c(6, 19)] <- c("beta.bmi", "se.bmi")


# merge smk, bmi, and t2d
temp_ces <- merge(select(smk_ces_clumped, SNP, beta.smk, se.smk), 
              select(pad_data, SNP, beta.outcome, se.outcome), by="SNP")
nrow(temp_ces)
combined_ces <- merge(select(temp_ces, SNP, beta.smk, se.smk, beta.outcome, se.outcome), 
                  select(bmi_ces_clumped, SNP, beta.bmi, se.bmi), by="SNP")
nrow(combined_ces)

# perform mvmr
mvmr_in_ces <- format_mvmr(BXGs = select(combined_ces, beta.smk, beta.bmi), 
                       BYG = combined_ces$beta.outcome, 
                       seBXGs = select(combined_ces, se.smk, se.bmi),
                       seBYG = combined_ces$se.outcome,
                       RSID = combined_ces$SNP)

# ******note: gencov is set to default 0. 
# ******may need  covariance between the effect of the genetic variants on each exposure
mvmr_out_ces <- mvmr(mvmr_in_ces, gencov = cov(combined_ces$beta.smk, combined_ces$beta.bmi), weights = 1)




#### try univariate MR using the same instrument: smoking initiation ####
uni <- smk_int[smk_int$SNP %in% combined$SNP,]

# LD clumping
uni_clumped <- clump_data(uni, clump_r2 = 0.1, clump_kb = 1)

# harmonising
uni_dat <- harmonise_data(exposure_dat = uni_clumped , outcome_dat = pad_data)

uni_res <- mr(uni_dat, method_list = c("mr_ivw", "mr_ivw_mre", "mr_ivw_fe", 
                                               "mr_egger_regression", "mr_egger_regression_bootstrap",
                                               "mr_weighted_median", "mr_simple_mode",
                                               "mr_weighted_mode"))

write.table(uni_res,"SMK.T2D.UNI.txt",quote=F,col.names=T,row.names=F,sep="\t")
#heterogeneity test
het = mr_heterogeneity(uni_dat)
write.table(het,"SMK.T2D.UNI.heterogeneity.txt",quote=F,col.names=T,row.names=F,sep="\t")
#horiz pleiotropy test
hp = mr_pleiotropy_test(uni_dat)
write.table(hp,"SMK.T2D.UNI.pleiotropy.txt",quote=F,col.names=T,row.names=F,sep="\t")


#Plots

#Scatterplot
pdf("SMK.T2D.UNI.Scatterplot.pdf")
p1 <- mr_scatter_plot(uni_res, uni_dat)
p1[[1]]
dev.off()

#Forest plot
pdf("SMK.T2D.UNI.Forestplot.pdf")
res_single <- mr_singlesnp(uni_dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]
dev.off()

pdf("SMK.T2D.UNI.OtherForestPlot.pdf")
res_single <- mr_singlesnp(uni_dat, all_method=c("mr_ivw", "mr_two_sample_ml"))
p3 <- mr_forest_plot(res_single)
p3[[1]]
dev.off()

#Leave One Out Plot
pdf("SMK.T2D.UNI.LeaveOneOutplot.pdf")
res_loo <- mr_leaveoneout(uni_dat)
p4 <- mr_leaveoneout_plot(res_loo)
p4[[1]]
dev.off()


#### try univariate MR using the same instrument: BMI ####
bmi <- read.table("/Users/dr/Desktop/MR-SMK.T2D/BMI-Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt", header=T)
bmi$phenotype = "BMI"
bmi <- format_data(bmi, type="exposure", snps=bmi$SNP, header=T, 
                   snp_col="SNP", phenotype_col = "phenotype",
                   beta_col="BETA", se_col="SE",
                   pval_col="P", samplesize_col ="N",
                   eaf_col="Freq_Tested_Allele_in_HRS",
                   effect_allele_col="Tested_Allele", 
                   other_allele_col="Other_Allele")
uni_bmi <- bmi[bmi$SNP %in% combined$SNP,]

# LD clumping
uni_bmi_clumped <- clump_data(uni_bmi, clump_r2 = 0.1, clump_kb = 1)

# harmonising
uni_bmi_dat <- harmonise_data(exposure_dat = uni_bmi_clumped , outcome_dat = pad_data)

uni_bmi_res <- mr(uni_bmi_dat, method_list = c("mr_ivw", "mr_ivw_mre", "mr_ivw_fe", 
                                       "mr_egger_regression", "mr_egger_regression_bootstrap",
                                       "mr_weighted_median"))

write.table(uni_bmi_res,"BMI.T2D.UNI.txt",quote=F,col.names=T,row.names=F,sep="\t")
#heterogeneity test
het = mr_heterogeneity(uni_bmi_dat)
write.table(het,"BMI.T2D.UNI.heterogeneity.txt",quote=F,col.names=T,row.names=F,sep="\t")
#horiz pleiotropy test
hp = mr_pleiotropy_test(uni_bmi_dat)
write.table(hp,"BMI.T2D.UNI.pleiotropy.txt",quote=F,col.names=T,row.names=F,sep="\t")


#Plots

#Scatterplot
pdf("BMI.T2D.UNI.Scatterplot.pdf")
p1 <- mr_scatter_plot(uni_bmi_res, uni_bmi_dat)
p1[[1]]
dev.off()

#Forest plot
pdf("BMI.T2D.UNI.Forestplot.pdf")
res_single <- mr_singlesnp(uni_bmi_dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]
dev.off()

pdf("BMI.T2D.UNI.OtherForestPlot.pdf")
res_single <- mr_singlesnp(uni_bmi_dat, all_method=c("mr_ivw", "mr_two_sample_ml"))
p3 <- mr_forest_plot(res_single)
p3[[1]]
dev.off()

#Leave One Out Plot
pdf("BMI.T2D.UNI.LeaveOneOutplot.pdf")
res_loo <- mr_leaveoneout(uni_bmi_dat)
p4 <- mr_leaveoneout_plot(res_loo)
p4[[1]]
dev.off()




#### univarite MR: smoking initiation - BMI ####
bmi <- read.table("/Users/dr/Desktop/MR-SMK.T2D/BMI-Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt", header=T)
bmi <- format_data(bmi, type="outcome", snps=bmi$SNP, header=T, 
                   snp_col="SNP", phenotype_col = "phenotype",
                   beta_col="BETA", se_col="SE",
                   pval_col="P", samplesize_col ="N",
                   eaf_col="Freq_Tested_Allele_in_HRS",
                   effect_allele_col="Tested_Allele", 
                   other_allele_col="Other_Allele")
uni_bmi <- bmi[bmi$SNP %in% combined$SNP,]

smk_int <- read_exposure_data("smoking_initiation/smoking_initiation.csv", sep=",")
smk_int_sub <- smk_int[smk_int$SNP %in% uni_bmi$SNP,]


# LD clumping
smk_int_sub_clumped <- clump_data(smk_int_sub, clump_r2 = 0.1, clump_kb = 1)

# harmonising
smk_int_sub_dat <- harmonise_data(exposure_dat = smk_int_sub_clumped, outcome_dat = uni_bmi)

smk_int_sub_res <- mr(smk_int_sub_dat, method_list = c("mr_ivw", "mr_ivw_mre", "mr_ivw_fe", 
                                               "mr_egger_regression", "mr_egger_regression_bootstrap",
                                               "mr_weighted_median"))

write.table(smk_int_sub_res,"SMK.BMI.txt",quote=F,col.names=T,row.names=F,sep="\t")
#heterogeneity test
het = mr_heterogeneity(smk_int_sub_dat)
write.table(het,"SMK.BMI.heterogeneity.txt",quote=F,col.names=T,row.names=F,sep="\t")
#horiz pleiotropy test
hp = mr_pleiotropy_test(smk_int_sub_dat)
write.table(hp,"SMK.BMI.pleiotropy.txt",quote=F,col.names=T,row.names=F,sep="\t")


#Plots

#Scatterplot
pdf("SMK.BMI.Scatterplot.pdf")
p1 <- mr_scatter_plot(smk_int_sub_res, smk_int_sub_dat)
p1[[1]]
dev.off()

#Forest plot
pdf("SMK.BMI.Forestplot.pdf")
res_single <- mr_singlesnp(smk_int_sub_dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]
dev.off()

pdf("SMK.BMI.OtherForestPlot.pdf")
res_single <- mr_singlesnp(smk_int_sub_dat, all_method=c("mr_ivw", "mr_two_sample_ml"))
p3 <- mr_forest_plot(res_single)
p3[[1]]
dev.off()

#Leave One Out Plot
pdf("SMK.BMI.LeaveOneOutplot.pdf")
res_loo <- mr_leaveoneout(smk_int_sub_dat)
p4 <- mr_leaveoneout_plot(res_loo)
p4[[1]]
dev.off()





#### different clumping ####
smk_r0.001 <- clump_data(smk_int, clump_r2 = 0.001, clump_kb = 1)
bmi_r0.001 <- clump_data(bmi_sub, clump_r2 = 0.001, clump_kb = 1)

smk_r0.1 <- clump_data(smk_int, clump_r2 = 0.1, clump_kb = 1)
bmi_r0.1 <- clump_data(bmi_sub, clump_r2 = 0.1, clump_kb = 1)


table(smk_r0.001 == smk_r0.1, useNA = 'ifany')
table(bmi_r0.001 == bmi_r0.1, useNA = 'ifany')






