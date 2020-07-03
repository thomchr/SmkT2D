


smk_clumped$r2 <- 2*(smk_clumped$beta.exposure)^2 * smk_clumped$eaf.exposure * (1-smk_clumped$eaf.exposure) /
  (2*(smk_clumped$beta.exposure)^2 * smk_clumped$eaf.exposure * (1-smk_clumped$eaf.exposure) + 
     (smk_clumped$se.exposure)^2*2*smk_clumped$samplesize.exposure*smk_clumped$eaf.exposure * (1-smk_clumped$eaf.exposure))
sum(smk_clumped$r2) # 0.009059311
