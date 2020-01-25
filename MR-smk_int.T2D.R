library(TwoSampleMR)
library(MRPRESSO)
library(ggplot2)

setwd("/Users/dr/Desktop/MR-SMK.T2D")

# load smk data
smk_int <-read_exposure_data("smoking_initiation/smoking_initiation.csv", sep=",")

# already loaded
# load t2d data
pad <- read.table("/Users/dr/Desktop/MR-SMK.T2D/Mahajan.NatGenet2018b.T2D.European.hg19rsid.txt", header=T)
pad$Phenotype <- "T2D"
pad_data <- format_data(pad, type="outcome", snps=pad$rsid, header=T, 
                        snp_col="rsid",
                        beta_col="Beta", se_col="SE",
                        eaf_col="EAF", effect_allele_col="EA",
                        other_allele_col="NEA", pval_col="Pvalue",
                        phenotype_col = "Phenotype")

# LD clumping
smk_int_clumped <- clump_data(smk_int, clump_r2 = 0.01, clump_kb = 250)


#### actual test SMK-T2D ####

#harmonize data
dat <- harmonise_data(exposure_dat = smk_int_clumped, outcome_dat = pad_data)
#run MR
res <- mr(dat, method_list = c("mr_ivw", "mr_ivw_mre", "mr_ivw_fe", 
                               "mr_egger_regression", "mr_egger_regression_bootstrap",
                               "mr_weighted_median", "mr_simple_mode",
                               "mr_weighted_mode"))

#res <- mr(dat, method_list = c("mr_ivw"))

#generate HTML file
#mr_report(dat)
#print output
write.table(res,"SMK.T2D.MR.clumped.txt",quote=F,col.names=T,row.names=F,sep="\t")
#heterogeneity test
het = mr_heterogeneity(dat)
write.table(het,"SMK.T2D.MR.clumped.heterogeneity.txt",quote=F,col.names=T,row.names=F,sep="\t")
#horiz pleiotropy test
hp = mr_pleiotropy_test(dat)
write.table(hp,"SMK.T2D.MR.clumped.pleiotropy.txt",quote=F,col.names=T,row.names=F,sep="\t")


#Plots

#Scatterplot
pdf("SMK.T2D.MR.clumped.Scatterplot.pdf")
p1 <- mr_scatter_plot(res, dat)
p1[[1]]
dev.off()

#Forest plot
pdf("SMK.T2D.MR.clumped.Forestplot.pdf")
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]
dev.off()

pdf("SMK.T2D.MR.clumped.OtherForestPlot.pdf")
res_single <- mr_singlesnp(dat, all_method=c("mr_ivw", "mr_two_sample_ml"))
p3 <- mr_forest_plot(res_single)
p3[[1]]
dev.off()

#Leave One Out Plot
pdf("SMK.T2D.MR.clumped.LeaveOneOutplot.pdf")
res_loo_ivw <- mr_leaveoneout(dat)
p4 <- mr_leaveoneout_plot(res_loo)
p4[[1]]
dev.off()

pdf("mr_egger_regression.SMK.T2D.MR.clumped.LeaveOneOutplot.pdf")
res_loo_egger <- mr_leaveoneout(dat, method = mr_egger_regression)
p_egger <- mr_leaveoneout_plot(res_loo)
p_egger[[1]]
dev.off()

#res_single <- mr_singlesnp(dat)
#p5 <- mr_funnel_plot(res_single)
#p5[[1]]
#dev.off()


#### MR-PRESSO ####
smk_int_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
          SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, 
          DISTORTIONtest = TRUE, data = dat, NbDistribution = 10000,  
          SignifThreshold = 0.05)

# remove outliers
smk_int_dist_outliers <- smk_int_presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
dat_outlier_removed <- dat[-smk_int_dist_outliers,]

res <- mr(dat_outlier_removed, method_list = c("mr_ivw", "mr_ivw_mre", "mr_ivw_fe", 
                                               "mr_egger_regression", "mr_egger_regression_bootstrap",
                                               "mr_weighted_median", "mr_simple_mode",
                                               "mr_weighted_mode"))
#res <- mr(dat_outlier_removed, method_list = c("mr_ivw"))

#generate HTML file
#mr_report(dat)
#print output
write.table(res,"presso.SMK.T2D.MR.clumped.txt",quote=F,col.names=T,row.names=F,sep="\t")
#heterogeneity test
het = mr_heterogeneity(dat_outlier_removed)
write.table(het,"presso.SMK.T2D.MR.clumped.heterogeneity.txt",quote=F,col.names=T,row.names=F,sep="\t")
#horiz pleiotropy test
hp = mr_pleiotropy_test(dat_outlier_removed)
write.table(hp,"presso.SMK.T2D.MR.clumped.pleiotropy.txt",quote=F,col.names=T,row.names=F,sep="\t")


#Plots

#Scatterplot
pdf("presso.SMK.T2D.MR.clumped.Scatterplot.pdf")
p1 <- mr_scatter_plot(res, dat_outlier_removed)
p1[[1]]
dev.off()

#Forest plot
pdf("presso.SMK.T2D.MR.clumped.Forestplot.pdf")
res_single <- mr_singlesnp(dat_outlier_removed)
p2 <- mr_forest_plot(res_single)
p2[[1]]
dev.off()

pdf("presso.SMK.T2D.MR.clumped.OtherForestPlot.pdf")
res_single <- mr_singlesnp(dat_outlier_removed, all_method=c("mr_ivw", "mr_two_sample_ml"))
p3 <- mr_forest_plot(res_single)
p3[[1]]
dev.off()

#Leave One Out Plot
pdf("ivw.presso.T2D.PAD.MR.clumped.LeaveOneOutplot.pdf")
res_loo_ivw <- mr_leaveoneout(dat_outlier_removed)
p4 <- mr_leaveoneout_plot(res_loo)
p4[[1]]
dev.off()

pdf("mr_egger_regression.presso.SMK.PAD.MR.clumped.LeaveOneOutplot.pdf")
res_loo_egger <- mr_leaveoneout(dat_outlier_removed, method = mr_egger_regression)
p_egger <- mr_leaveoneout_plot(res_loo)
p_egger[[1]]
dev.off()

#res_single <- mr_singlesnp(dat)
#p5 <- mr_funnel_plot(res_single)
#p5[[1]]
#dev.off()




#### positive control SMK-CAD ####

#cad <- read.table("/Users/dr/Desktop/MR-SMK.T2D/CAD_META", header=T)
#cad$Phenotype <- "CAD"
#cad_data <- format_data(cad, type="outcome", snps=cad$oldID, header=T, 
#                        snp_col="oldID",
#                        beta_col="Effect", se_col="StdErr",
#                        eaf_col="Feq1", effect_allele_col="Allele1",
#                        other_allele_col="Allele2", pval_col="P-value",
#                        phenotype_col = "Phenotype")


#harmonize data
cad_dat <- harmonise_data(exposure_dat = smk_int_clumped, outcome_dat = cad_data)
#run MR
res <- mr(cad_dat, method_list = c("mr_ivw", "mr_ivw_mre", "mr_ivw_fe", 
                               "mr_egger_regression", "mr_egger_regression_bootstrap",
                               "mr_weighted_median", "mr_simple_mode",
                               "mr_weighted_mode"))
#res <- mr(cad_dat, method_list = c("mr_ivw"))
#generate HTML file
#mr_report(dat)
#print output
write.table(res,"SMK.CAD.MR.clumped.txt",quote=F,col.names=T,row.names=F,sep="\t")
#heterogeneity test
het = mr_heterogeneity(cad_dat)
write.table(het,"SMK.CAD.MR.clumped.heterogeneity.txt",quote=F,col.names=T,row.names=F,sep="\t")
#horiz pleiotropy test
hp = mr_pleiotropy_test(cad_dat)
write.table(hp,"SMK.CAD.MR.clumped.pleiotropy.txt",quote=F,col.names=T,row.names=F,sep="\t")


#Plots

#Scatterplot
pdf("SMK.CAD.MR.clumped.Scatterplot.pdf")
p1 <- mr_scatter_plot(res, cad_dat)
p1[[1]]
dev.off()

#Forest plot
pdf("SMK.CAD.MR.clumped.Forestplot.pdf")
res_single <- mr_singlesnp(cad_dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]
dev.off()

pdf("SMK.CAD.MR.clumped.OtherForestPlot.pdf")
res_single <- mr_singlesnp(dat, all_method=c("mr_ivw", "mr_two_sample_ml"))
p3 <- mr_forest_plot(res_single)
p3[[1]]
dev.off()

#Leave One Out Plot
pdf("SMK.CAD.MR.clumped.LeaveOneOutplot.pdf")
res_loo <- mr_leaveoneout(cad_dat)
p4 <- mr_leaveoneout_plot(res_loo)
p4[[1]]
dev.off()

#res_single <- mr_singlesnp(dat)
#p5 <- mr_funnel_plot(res_single)
#p5[[1]]
#dev.off()

#### positive control CAD MR-presso ####

#cad <- read.table("/Users/dr/Desktop/MR-SMK.T2D/CAD_META", header=T)
#cad_data <- format_data(cad, type="outcome", snps=cad$oldID, header=T, 
#                        snp_col="oldID",
#                        beta_col="Effect", se_col="StdErr",
#                        eaf_col="Feq1", effect_allele_col="Allele1",
#                        other_allele_col="Allele2", pval_col="P-value")


#harmonize data
dat <- harmonise_data(exposure_dat = smk_int_clumped, outcome_dat = cad_data)

smk_int_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                            SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, 
                            DISTORTIONtest = TRUE, data = dat, NbDistribution = 10000,  
                            SignifThreshold = 0.05)

# remove outliers
smk_int_dist_outliers <- smk_int_presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
dat_outlier_removed <- dat[-smk_int_dist_outliers,]

#run MR
#res <- mr(dat)
res <- mr(dat_outlier_removed, method_list = c("mr_ivw", "mr_ivw_mre", "mr_ivw_fe", 
                                               "mr_egger_regression", "mr_egger_regression_bootstrap",
                                               "mr_weighted_median", "mr_simple_mode",
                                               "mr_weighted_mode"))


#generate HTML file
#mr_report(dat)
#print output
write.table(res,"presso.SMK.CAD.MR.clumped.txt",quote=F,col.names=T,row.names=F,sep="\t")
#heterogeneity test
het = mr_heterogeneity(dat_outlier_removed)
write.table(het,"presso.SMK.CAD.MR.clumped.heterogeneity.txt",quote=F,col.names=T,row.names=F,sep="\t")
#horiz pleiotropy test
hp = mr_pleiotropy_test(dat_outlier_removed)
write.table(hp,"presso.SMK.CAD.MR.clumped.pleiotropy.txt",quote=F,col.names=T,row.names=F,sep="\t")


#Plots

#Scatterplot
pdf("presso.SMK.CAD.MR.clumped.Scatterplot.pdf")
p1 <- mr_scatter_plot(res, dat_outlier_removed)
p1[[1]]
dev.off()

#Forest plot
pdf("presso.SMK.CAD.MR.clumped.Forestplot.pdf")
res_single <- mr_singlesnp(dat_outlier_removed)
p2 <- mr_forest_plot(res_single)
p2[[1]]
dev.off()

pdf("presso.SMK.CAD.MR.clumped.OtherForestPlot.pdf")
res_single <- mr_singlesnp(dat_outlier_removed, all_method=c("mr_ivw", "mr_two_sample_ml"))
p3 <- mr_forest_plot(res_single)
p3[[1]]
dev.off()

#Leave One Out Plot
pdf("presso.SMK.CAD.MR.clumped.LeaveOneOutplot.pdf")
res_loo <- mr_leaveoneout(dat_outlier_removed)
p4 <- mr_leaveoneout_plot(res_loo)
p4[[1]]
dev.off()

#res_single <- mr_singlesnp(dat)
#p5 <- mr_funnel_plot(res_single)
#p5[[1]]
#dev.off()



