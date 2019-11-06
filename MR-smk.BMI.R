library(MVMR)
library(TwoSampleMR)
library(dplyr)
setwd("/Users/dr/Desktop/MR-SMK.T2D")

#### BMI -> SMK ####


# load smoking data as 
smk_int <-read_outcome_data("smoking_initiation/smoking_initiation.csv", sep=",")

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
#bmi_sub <- bmi[bmi$SNP %in% combined$SNP,]

# clumping
bmi_clumped <- clump_data(bmi_sub, clump_r2 = 0.1, clump_kb = 1)

# harmonising
bmi_harm <- harmonise_data(exposure_dat = bmi_clumped, outcome_dat = smk_int)
nrow(bmi_harm)

bmi_smk_res <- mr(bmi_harm, method_list = c("mr_ivw", "mr_ivw_mre", "mr_ivw_fe", 
                                                       "mr_egger_regression", "mr_egger_regression_bootstrap",
                                                       "mr_weighted_median"))

write.table(bmi_smk_res,"BMI.SMK.txt",quote=F,col.names=T,row.names=F,sep="\t")
#heterogeneity test
het = mr_heterogeneity(bmi_harm)
write.table(het,"BMI.SMK.heterogeneity.txt",quote=F,col.names=T,row.names=F,sep="\t")
#horiz pleiotropy test
hp = mr_pleiotropy_test(bmi_harm)
write.table(hp,"BMI.SMK.pleiotropy.txt",quote=F,col.names=T,row.names=F,sep="\t")


#Plots

#Scatterplot
pdf("BMI.SMK.Scatterplot.pdf")
p1 <- mr_scatter_plot(bmi_smk_res, bmi_harm)
p1[[1]]
dev.off()

#Forest plot
pdf("BMI.SMK.Forestplot.pdf")
res_single <- mr_singlesnp(bmi_harm)
p2 <- mr_forest_plot(res_single)
p2[[1]]
dev.off()

pdf("BMI.SMK.OtherForestPlot.pdf")
res_single <- mr_singlesnp(bmi_harm, all_method=c("mr_ivw", "mr_two_sample_ml"))
p3 <- mr_forest_plot(res_single)
p3[[1]]
dev.off()

#Leave One Out Plot
pdf("BMI.SMK.LeaveOneOutplot.pdf")
res_loo <- mr_leaveoneout(bmi_harm)
p4 <- mr_leaveoneout_plot(res_loo)
p4[[1]]
dev.off()







#### SMK -> BMI ####
# load smoking data as 
smk_int <-read_exposure_data("smoking_initiation/smoking_initiation.csv", sep=",")

# bmi data
bmi <- read.table("/Users/dr/Desktop/MR-SMK.T2D/BMI-Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt", header=T)
bmi$phenotype = "BMI"
bmi <- format_data(bmi, type="outcome", snps=bmi$SNP, header=T, 
                   snp_col="SNP", phenotype_col = "phenotype",
                   beta_col="BETA", se_col="SE",
                   pval_col="P", samplesize_col ="N",
                   eaf_col="Freq_Tested_Allele_in_HRS",
                   effect_allele_col="Tested_Allele", 
                   other_allele_col="Other_Allele")
smk_int_sub <- smk_int[smk_int$SNP %in% combined$SNP,]

# clumping
smk_clumped <- clump_data(smk_int_sub, clump_r2 = 0.01, clump_kb = 250)

# harmonising
smk_harm <- harmonise_data(exposure_dat = smk_clumped, outcome_dat = bmi)
nrow(smk_harm)

smk_bmi_res <- mr(smk_harm, method_list = c("mr_ivw", "mr_ivw_mre", "mr_ivw_fe", 
                                            "mr_egger_regression", "mr_egger_regression_bootstrap",
                                            "mr_weighted_median"))

write.table(smk_bmi_res,"SMK.BMI.txt",quote=F,col.names=T,row.names=F,sep="\t")
#heterogeneity test
het = mr_heterogeneity(smk_harm)
write.table(het,"SMK.BMI.heterogeneity.txt",quote=F,col.names=T,row.names=F,sep="\t")
#horiz pleiotropy test
hp = mr_pleiotropy_test(smk_harm)
write.table(hp,"SMK.BMI.pleiotropy.txt",quote=F,col.names=T,row.names=F,sep="\t")


#Plots

#Scatterplot
pdf("SMK.BMI.Scatterplot.pdf")
p1 <- mr_scatter_plot(smk_bmi_res, smk_harm)
p1[[1]]
dev.off()

#Forest plot
pdf("SMK.BMI.Forestplot.pdf")
res_single <- mr_singlesnp(smk_harm)
p2 <- mr_forest_plot(res_single)
p2[[1]]
dev.off()

pdf("SMK.BMI.OtherForestPlot.pdf")
res_single <- mr_singlesnp(smk_harm, all_method=c("mr_ivw", "mr_two_sample_ml"))
p3 <- mr_forest_plot(res_single)
p3[[1]]
dev.off()

#Leave One Out Plot
pdf("SMK.BMI.LeaveOneOutplot.pdf")
res_loo <- mr_leaveoneout(smk_harm)
p4 <- mr_leaveoneout_plot(res_loo)
p4[[1]]
dev.off()





