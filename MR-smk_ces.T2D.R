library(TwoSampleMR)

setwd("/Users/dr/Desktop/MR-SMK.T2D")

# load smk data
smk_ces <-read_exposure_data("smoking_cessation/smk_ces.csv", sep=",")

# already loaded
# load t2d data
#pad <- read.table("/Users/dr/Desktop/MR-SMK.T2D/Mahajan.NatGenet2018b.T2D.European.hg19rsid.txt", header=T)
#pad_data <- format_data(pad, type="outcome", snps=pad$rsid, header=T, 
#                        snp_col="rsid",
#                        beta_col="Beta", se_col="SE",
#                        eaf_col="EAF", effect_allele_col="EA",
#                        other_allele_col="NEA", pval_col="Pvalue",
#                        phenotype_col = "Phenotype")

# LD clumping
smk_ces_clumped <- clump_data(smk_ces, clump_r2 = 0.01, clump_kb = 250)


#### actual test SMK-T2D ####

#harmonize data
smk_ces_dat <- harmonise_data(exposure_dat = smk_ces_clumped, outcome_dat = pad_data)
#run MR
smk_ces_res <- mr(smk_ces_dat, method_list = c("mr_ivw", "mr_ivw_mre", "mr_ivw_fe", 
                                               "mr_egger_regression", "mr_egger_regression_bootstrap",
                                               "mr_weighted_median", "mr_simple_mode",
                                               "mr_weighted_mode"))
#smk_ces_res <- mr(smk_ces_dat, method_list = c("mr_ivw"))



#generate HTML file
#mr_report(dat)
#print output
write.table(smk_ces_res,"SMK.T2D.MR.clumped.txt",quote=F,col.names=T,row.names=F,sep="\t")
#heterogeneity test
het = mr_heterogeneity(smk_ces_dat)
write.table(het,"SMK.T2D.MR.clumped.heterogeneity.txt",quote=F,col.names=T,row.names=F,sep="\t")
#horiz pleiotropy test
hp = mr_pleiotropy_test(smk_ces_dat)
write.table(hp,"SMK.T2D.MR.clumped.pleiotropy.txt",quote=F,col.names=T,row.names=F,sep="\t")


#Plots

#Scatterplot
pdf("SMK.T2D.MR.clumped.Scatterplot.pdf")
p1 <- mr_scatter_plot(smk_ces_res, smk_ces_dat)
p1[[1]]
dev.off()

#Forest plot
pdf("SMK.T2D.MR.clumped.Forestplot.pdf")
res_single <- mr_singlesnp(smk_ces_dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]
dev.off()

pdf("SMK.T2D.MR.clumped.OtherForestPlot.pdf")
res_single <- mr_singlesnp(smk_ces_dat, all_method=c("mr_ivw", "mr_two_sample_ml"))
p3 <- mr_forest_plot(res_single)
p3[[1]]
dev.off()

#Leave One Out Plot
pdf("SMK.T2D.MR.clumped.LeaveOneOutplot.pdf")
res_loo <- mr_leaveoneout(smk_ces_dat)
p4 <- mr_leaveoneout_plot(res_loo)
p4[[1]]
dev.off()

#res_single <- mr_singlesnp(dat)
#p5 <- mr_funnel_plot(res_single)
#p5[[1]]
#dev.off()



#### CAD ####
#cad <- read.table("/Users/dr/Desktop/MR-SMK.T2D/CAD_META", header=T)
#cad$Phenotype <- "CAD"
#cad_data <- format_data(cad, type="outcome", snps=cad$oldID, header=T, 
#                        snp_col="oldID",
#                        beta_col="Effect", se_col="StdErr",
#                        eaf_col="Feq1", effect_allele_col="Allele1",
#                        other_allele_col="Allele2", pval_col="P-value",
#                        phenotype_col = "Phenotype")



#harmonize data
cad_smk_ces_dat <- harmonise_data(exposure_dat = smk_ces_clumped, outcome_dat = cad_data)
#run MR
#cad_smk_ces_res <- mr(cad_smk_ces_dat, method_list = c("mr_ivw", "mr_ivw_mre", "mr_ivw_fe", 
#                               "mr_egger_regression", "mr_egger_regression_bootstrap",
#                               "mr_weighted_median", "mr_simple_mode",
#                               "mr_weighted_mode"))
cad_smk_ces_res <- mr(cad_smk_ces_dat, method_list = c("mr_ivw"))

#generate HTML file
#mr_report(dat)
#print output
write.table(cad_smk_ces_res,"SMK.CAD.MR.clumped.txt",quote=F,col.names=T,row.names=F,sep="\t")
#heterogeneity test
het = mr_heterogeneity(cad_smk_ces_dat)
write.table(het,"SMK.CAD.MR.clumped.heterogeneity.txt",quote=F,col.names=T,row.names=F,sep="\t")
#horiz pleiotropy test
hp = mr_pleiotropy_test(cad_smk_ces_dat)
write.table(hp,"SMK.CAD.MR.clumped.pleiotropy.txt",quote=F,col.names=T,row.names=F,sep="\t")


#Plots

#Scatterplot
pdf("SMK.CAD.MR.clumped.Scatterplot.pdf")
p1 <- mr_scatter_plot(cad_smk_ces_res, cad_smk_ces_dat)
p1[[1]]
dev.off()

#Forest plot
pdf("SMK.CAD.MR.clumped.Forestplot.pdf")
res_single <- mr_singlesnp(cad_smk_ces_dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]
dev.off()

pdf("SMK.CAD.MR.clumped.OtherForestPlot.pdf")
res_single <- mr_singlesnp(cad_smk_ces_dat, all_method=c("mr_ivw", "mr_two_sample_ml"))
p3 <- mr_forest_plot(res_single)
p3[[1]]
dev.off()

#Leave One Out Plot
pdf("SMK.CAD.MR.clumped.LeaveOneOutplot.pdf")
res_loo <- mr_leaveoneout(cad_smk_ces_dat)
p4 <- mr_leaveoneout_plot(res_loo)
p4[[1]]
dev.off()

#res_single <- mr_singlesnp(dat)
#p5 <- mr_funnel_plot(res_single)
#p5[[1]]
#dev.off()

