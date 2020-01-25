#### smoking initiation - T2D ####
m.smk_int <- merge(x = dat, y = pad, by.x = "SNP", by.y = "rsid")


m.smk_int.sub <- m.smk_int[, c("SNP", "Chr", "Pos", "effect_allele.exposure", "other_allele.exposure",
                               "eaf.exposure", "beta.exposure", "se.exposure", "beta.outcome", "se.outcome")]

colnames(m.smk_int.sub) <- c("rsid", "chromosome", "position", "effect allele", "non-effect allele",
                             "effect allele frequency", "smoking-initiation effect size", 
                             "smoking-initiation standard error", "T2D effect size", "T2D standard error")

write.csv(m.smk_int.sub, "/Users/dr/Desktop/MR.SMK-T2D/instrument.smoking_initiation-T2D.csv")


#### smoking initiation - CAD ####
m.smk_int <- merge(x = cad_dat, y = cad, by.x = "SNP", by.y = "oldID")


m.smk_int.sub <- m.smk_int[, c("SNP", "CHR", "BP", "effect_allele.exposure", "other_allele.exposure",
                               "eaf.exposure", "beta.exposure", "se.exposure", "beta.outcome", "se.outcome")]

colnames(m.smk_int.sub) <- c("rsid", "chromosome", "position", "effect allele", "non-effect allele",
                             "effect allele frequency", "smoking-initiation effect size", 
                             "smoking-initiation standard error", "T2D effect size", "T2D standard error")

write.csv(m.smk_int.sub, "/Users/dr/Desktop/MR.SMK-T2D/instrument.smoking_initiation-CAD.csv")

#### smoking cessation - T2D ####
m.smk_ces <- merge(x = smk_ces_dat, y = pad, by.x = "SNP", by.y = "rsid")


m.smk_ces.sub <- m.smk_ces[, c("SNP", "Chr", "Pos", "effect_allele.exposure", "other_allele.exposure",
                               "eaf.exposure", "beta.exposure", "se.exposure", "beta.outcome", "se.outcome")]

colnames(m.smk_ces.sub) <- c("rsid", "chromosome", "position", "effect allele", "non-effect allele",
                             "effect allele frequency", "smoking-cessation effect size", 
                             "smoking-cessation standard error", "T2D effect size", "T2D standard error")
write.csv(m.smk_ces.sub, "/Users/dr/Desktop/MR.SMK-T2D/instrument.smoking_cessation.csv")

#### cigarettes per day - T2D ####
m.cig_day <- merge(x = cig_day_dat, y = pad, by.x = "SNP", by.y = "rsid")


m.cig_day.sub <- m.cig_day[, c("SNP", "Chr", "Pos", "effect_allele.exposure", "other_allele.exposure",
                               "eaf.exposure", "beta.exposure", "se.exposure", "beta.outcome", "se.outcome")]

colnames(m.cig_day.sub) <- c("rsid", "chromosome", "position", "effect allele", "non-effect allele",
                             "effect allele frequency", "cigarettes-per-day effect size", 
                             "cigarettes-per-day standard error", "T2D effect size", "T2D standard error")
write.csv(m.cig_day.sub, "/Users/dr/Desktop/MR.SMK-T2D/instrument.cigarettes-per-day.csv")
