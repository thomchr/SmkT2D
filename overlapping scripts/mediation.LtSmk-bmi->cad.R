
mvmr.cad <- read.csv("/Users/dr/Desktop/MR.SMK-T2D/better_instrumnet_mvmr/instruments/combined_IV.LtSmk_bmi->cad.instrument.csv", 
                     stringsAsFactors = FALSE)


# setups
# assume idenpendent SNPs after clumping
nrow(mvmr.cad) # 1441
rho <- diag(nrow(mvmr.cad))

# snp -> risk factor
betaXG <- mvmr.cad$lifetime.smoking.effect.size
sebetaXG <- mvmr.cad$lifetime.smoking.standard.error

# snp -> mediator
betaMG <- mvmr.cad$BMI.effect.size
sebetaMG <- mvmr.cad$BMI.standard.error

# snp -> outcome
betaYG <- mvmr.cad$CAD.effect.size
sebetaYG <- mvmr.cad$CAD.standard.error


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


result.cad <- data.frame("total_effect" = total.effect.correl,
                         "se.total.fixed" = se.total.effect.fixed,
                         "se.total.random" = se.total.effect.random,
                         "direct_effect" = direct.effect.correl,
                         "se.direct.fixed" = se.direct.effect.fixed,
                         "se.direct.random" = se.direct.effect.random)
write.csv(result.cad, "/Users/dr/Desktop/MR.SMK-T2D/better_instrumnet_mvmr/combined_IV.direct_total.cad.csv")  







