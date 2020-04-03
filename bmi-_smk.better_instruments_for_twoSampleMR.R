#### bmi instrument ####

#### overlaps ####
# load smoking data
smk_int <- read_exposure_data("smoking_initiation/smoking_initiation.csv", sep=",")

# load t2d data
pad <- read.table("/Users/dr/Desktop/MR.SMK-T2D/Mahajan.NatGenet2018b.T2D.European.hg19rsid.txt", header=T)
pad_data <- format_data(pad, type="outcome", snps=pad$rsid, header=T,
                        snp_col="rsid",
                        beta_col="Beta", se_col="SE",
                        eaf_col="EAF", effect_allele_col="EA",
                        other_allele_col="NEA", pval_col="Pvalue")


# load bmi data
bmi <- read.table("/Users/dr/Desktop/MR.SMK-T2D/BMI-Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt", header=T)
bmi$phenotype = "BMI"


# only keep genome-wide significant snps
bmi <- bmi[bmi$P <= 5*10^-8,]
bmi <- format_data(bmi, type="exposure", snps=bmi$SNP, header=T,
                   snp_col="SNP", phenotype_col = "phenotype",
                   beta_col="BETA", se_col="SE",
                   pval_col="P", samplesize_col ="N",
                   eaf_col="Freq_Tested_Allele_in_HRS",
                   effect_allele_col="Tested_Allele",
                   other_allele_col="Other_Allele",
                   chr_col="CHR", pos_col="POS")

# clumping bmi 
bmi_sub_clumped <- clump_data(bmi, clump_r2 = 0.01, clump_kb = 250)

# harmonising by T2D data
smk_int_harm <- harmonise_data(exposure_dat = smk_int, outcome_dat = pad_data)
bmi_harm <- harmonise_data(exposure_dat = bmi_sub_clumped, outcome_dat = pad_data)

# merge smk and bmi instruments
combined <- merge(select(smk_int_harm, SNP, effect_allele.exposure, other_allele.exposure, eaf.exposure,
                         beta.exposure, se.exposure, pval.exposure), 
                  select(bmi_harm, SNP, effect_allele.exposure, other_allele.exposure, eaf.exposure,
                         beta.exposure, se.exposure, pval.exposure, chr.exposure, pos.exposure, samplesize.exposure), by="SNP")

#### load smk proxies ####
# load smk proxies
bmi_smk.px <- read.table("/Users/dr/Desktop/MR.SMK-T2D/BMI_snp_proxies/Proxies.BMISmk.txt", header = TRUE)

# seperate out smk for harmonising
smk.px <- data.frame("SNP" = bmi_smk.px$SNP_B, "effect_allele.exposure" = bmi_smk.px$ALT,
                     "other_allele.exposure" = bmi_smk.px$REF, "eaf.exposure" = bmi_smk.px$AF,
                     "beta.exposure" = bmi_smk.px$BETA.1, "se.exposure" = bmi_smk.px$SE.1,
                     "pval.exposure" = bmi_smk.px$PVALUE, "samplesize.exposure" = bmi_smk.px$N.1,
                     "exposure" = "smoking", "mr_keep.exposure" = TRUE, 
                     "pval_origin.exposure" = "reported", "id.exposure" = "qYo4wL",
                     "name" = bmi_smk.px$SNP_B)

# seperate out bmi for harmonising
bmi.px <- data.frame("SNP" = bmi_smk.px$SNP_B, "effect_allele.exposure" = bmi_smk.px$Tested_Allele,
                     "other_allele.exposure" = bmi_smk.px$Other_Allele, 
                     "eaf.exposure" = bmi_smk.px$Freq_Tested_Allele_in_HRS,
                     "beta.exposure" = bmi_smk.px$BETA, "se.exposure" = bmi_smk.px$SE,
                     "pval.exposure" = bmi_smk.px$P, "samplesize.exposure" = bmi_smk.px$N,
                     "exposure" = "bmi", "mr_keep.exposure" = TRUE, 
                     "pval_origin.exposure" = "reported", "id.exposure" = "uLruT0",
                     "name" = bmi_smk.px$SNP_B,
                     "chr.exposure" = bmi_smk.px$CHR_B, "pos.exposure" = bmi_smk.px$BP_B)

# clumping
bmi.px.clumped <- clump_data(bmi.px, clump_r2 = 0.01, clump_kb = 250)
smk.px.clumped <- clump_data(smk.px, clump_r2 = 0.01, clump_kb = 250)

# harmonising by T2D data
bmi.px.harm <- harmonise_data(exposure_dat = bmi.px.clumped, outcome_dat = pad_data)
smk.px.harm <- harmonise_data(exposure_dat = smk.px.clumped, outcome_dat = pad_data)


# merge smk and bmi instruments
combined.px <- merge(select(smk.px.harm, SNP, beta.exposure, se.exposure, pval.exposure, 
                            effect_allele.exposure, other_allele.exposure, eaf.exposure), 
                     select(bmi.px.harm, SNP, beta.exposure, se.exposure, pval.exposure, 
                            effect_allele.exposure, other_allele.exposure, eaf.exposure,
                            chr.exposure, pos.exposure, samplesize.exposure), 
                     by="SNP")

BMISmk_instrument <- data.frame("SNP" = c(as.character(combined$SNP), as.character(combined.px$SNP)), 
                                "chr" = c(combined$chr.exposure, combined.px$chr.exposure),
                                "pos" = c(combined$pos.exposure, combined.px$pos.exposure),
                                "samplesize" = c(combined$samplesize.exposure, combined.px$samplesize.exposure),
                                
                                # smoking columns
                                "smk.effect_allele" = c(as.character(combined$effect_allele.exposure.x), 
                                                        as.character(combined.px$effect_allele.exposure.x)),
                                "smk.other_allele" = c(as.character(combined$other_allele.exposure.x), 
                                                       as.character(combined.px$other_allele.exposure.x)),
                                "smk.eaf" = c(combined$eaf.exposure.x, combined.px$eaf.exposure.x),
                                "smk.beta" = c(combined$beta.exposure.x, combined.px$beta.exposure.x),
                                "smk.se" = c(combined$se.exposure.x, combined.px$se.exposure.x),
                                "smk.pval" = c(combined$pval.exposure.x, combined.px$pval.exposure.x),
                                
                                # bmi columns
                                "bmi.effect_allele" = c(as.character(combined$effect_allele.exposure.y), 
                                                        as.character(combined.px$effect_allele.exposure.y)),
                                "bmi.other_allele" = c(as.character(combined$other_allele.exposure.y), 
                                                       as.character(combined.px$other_allele.exposure.y)),
                                "bmi.eaf" = c(combined$eaf.exposure.y, combined.px$eaf.exposure.y),
                                "bmi.beta" = c(combined$beta.exposure.y, combined.px$beta.exposure.y),
                                "bmi.se" = c(combined$se.exposure.y, combined.px$se.exposure.y),
                                "bmi.pval" = c(combined$pval.exposure.y, combined.px$pval.exposure.y))
nrow(BMISmk_instrument)

write.csv(BMISmk_instrument, 
          "/Users/dr/Desktop/MR.SMK-T2D/better_instrument_twoSampleMR/BMISmk.instrument.csv",
          row.names = FALSE)












