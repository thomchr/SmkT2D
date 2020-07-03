
mvmr.t2d <- read.csv("/Users/dr/Desktop/MR.SMK-T2D/non_overlap_analyses/instruments/instrument.SmkInt_bmi->t2d.csv", 
                     stringsAsFactors = FALSE)


# setups
# assume idenpendent SNPs after clumping
nrow(mvmr.t2d) # 207
rho <- diag(nrow(mvmr.t2d))

# snp -> risk factor
betaXG <- mvmr.t2d$smoking.initiation.effect.size
sebetaXG <- mvmr.t2d$smoking.initiation.standard.error

# snp -> mediator
betaMG <- mvmr.t2d$BMI.effect.size
sebetaMG <- mvmr.t2d$BMI.standard.error

# snp -> outcome
betaYG <- mvmr.t2d$T2D.effect.size
sebetaYG <- mvmr.t2d$T2D.standard.error


# calculation
Omega <- sebetaYG%o%sebetaYG*rho

total.effect.correl <- solve(t(betaXG)%*%solve(Omega)%*%betaXG)*t(betaXG)%*%solve(Omega)%*%betaYG

se.total.effect.fixed <- sqrt(solve(t(betaXG)%*%solve(Omega)%*%betaXG))

resid.total <- betaYG-total.effect.correl*betaXG

se.total.effect.random <- sqrt(solve(t(betaXG)%*%solve(Omega)%*%betaXG))* 
  max(sqrt(t(resid.total)%*%solve(Omega)%*%resid.total/(length(betaXG)-1)),1)
direct.effect.correl <- solve(t(cbind(betaXG, betaMG))%*%solve(Omega)%*%
                                cbind(betaXG, betaMG))%*%t(cbind(betaXG, betaMG))%*%solve(Omega)%*%betaYG

se.direct.effect.fixed <- c(sqrt(solve(t(cbind(betaXG, betaMG))%*%solve(Omega)%*%cbind(betaXG, betaMG))[1,1]),
                            sqrt(solve(t(cbind(betaXG, betaMG))%*%solve(Omega)%*%cbind(betaXG, betaMG))[2,2]))

resid.direct <- betaYG-direct.effect.correl[1]*betaXG-direct.effect.correl[2]*betaMG

se.direct.effect.random <- c(sqrt(solve(t(cbind(betaXG, betaMG))%*%solve(Omega)%*%cbind(betaXG, betaMG))[1,1]), 
                             sqrt(solve(t(cbind(betaXG, betaMG))%*%solve(Omega)%*%cbind(betaXG, betaMG))[2,2]))*
  max(sqrt(t(resid.direct)%*%solve(Omega)%*%resid.direct/(length(betaXG)-2)),1)


result.t2d <- data.frame("total_effect" = total.effect.correl,
                         "se.total.fixed" = se.total.effect.fixed,
                         "se.total.random" = se.total.effect.random,
                         "direct_effect" = direct.effect.correl,
                         "se.direct.fixed" = se.direct.effect.fixed,
                         "se.direct.random" = se.direct.effect.random)
#write.csv(result.t2d, "/Users/dr/Desktop/MR.SMK-T2D/non_overlap_analyses/nonoverlap.combined_IV.mediation.t2d.csv")  







