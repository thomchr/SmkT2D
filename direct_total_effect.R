mvmr.t2d <- read.csv("/Users/dr/Desktop/MR.SMK-T2D/better_instrumnet_mvmr/instruments/mvmr_instrument.smk-bmi->t2d.csv", 
                     stringsAsFactors = FALSE)
mvmr.cad <- read.csv("/Users/dr/Desktop/MR.SMK-T2D/better_instrumnet_mvmr/instruments/mvmr_instrument.smk-bmi->cad.csv", 
                     stringsAsFactors = FALSE)
which(mvmr.t2d$RSID %in% mvmr.cad$RSID == FALSE)  # exactly same instrument

# setups
# assume idenpendent SNPs after clumping
rho <- diag(286)

# snp -> risk factor
betaXG <- mvmr.t2d$beta.smk
sebetaXG <- mvmr.t2d$se.smk

# snp -> mediator
betaMG <- mvmr.t2d$beta.bmi
sebetaMG <- mvmr.t2d$se.bmi

# snp -> outcome
betaYG <- mvmr.t2d$beta.t2d
sebetaYG <- mvmr.t2d$se.t2d


# calculation
Omega <- sebetaYG%o%sebetaYG*rho

total.effect.correl <- solve(t(betaXG)%*%solve(Omega)%*%betaXG)*t(betaXG)%*%solve(Omega)%*%betaYG

se.total.effect.fixed <- sqrt(solve(t(betaXG)%*%solve(Omega)%*%betaXG))

resid.total <- betaYG-total.effect.correl*betaXG

se.total.effect.random <- sqrt(solve(t(betaXG)%*%solve(Omega)%*%betaXG))* 
  max(sqrt(t(resid.total)%*%solve(Omega)%*%resid.total/(length(betaXG)-1)),1)
direct.effect.correl <- solve(t(cbind(betaXG, betaMG))%*%solve(Omega)%*%
                                cbind(betaXG, betaMG))%*%t(cbind(betaXG, betaMG))%*%solve(Omega)%*%betaYG
                                
se.direct.effect.fixed <- sqrt(solve(t(cbind(betaXG, betaMG))%*%solve(Omega)%*%cbind(betaXG, betaMG))[1,1])

resid.direct <- betaYG-direct.effect.correl[1]*betaXG-direct.effect.correl[2]*betaMG

se.direct.effect.random <- sqrt(solve(t(cbind(betaXG, betaMG))%*%solve(Omega)%*%cbind(betaXG, betaMG))[1,1])*
  max(sqrt(t(resid.direct)%*%solve(Omega)%*%resid.direct/(length(betaXG)-2)),1)
  
  
  





