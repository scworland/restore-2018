library(lmomco)
library(feather)
looPPLOmodel <- read_feather("../../../results/gage/gam/all_gage_looest_pplo.feather")
looL1model   <- read_feather("../../../results/gage/gam/all_gage_looest_L1.feather"  )
looT2model   <- read_feather("../../../results/gage/gam/all_gage_looest_T2.feather"  )
looT3model   <- read_feather("../../../results/gage/gam/all_gage_looest_T3.feather"  )
looT4model   <- read_feather("../../../results/gage/gam/all_gage_looest_T4.feather"  )
looT5model   <- read_feather("../../../results/gage/gam/all_gage_looest_T5.feather"  )




#ObsMean <- 0*looPPLOmodel$loo_est_pplo + (1-looPPLOmodel$loo_est_pplo)*looL1model$loo_est_L1
EstMeanFlow_loo <- (1-looPPLOmodel$loo_est_pplo)*looL1model$loo_est_L1*looL1model$loo_bias_corr
EstMeanFlow <- (1-looPPLOmodel$est_pplo)*looL1model$est_L1*looL1model$bias_corr
EstMeanFlow_loo[looL1model$site_no == "08167000"]
EstMeanFlow[looL1model$site_no == "08167000"]
plot(EstMeanFlow, EstMeanFlow_loo, log='xy', pch=16, col=rgb(1,0,0,.3), lwd=0.5, cex=1,
     xlab="Estimated Decadal Mean Flow", ylab="LOO Estimated Decadle Mean Flow")
abline(0,1)

ObsMean <- (1-looPPLOmodel$pplo)*looL1model$L1

median(abs(log10(EstMeanFlow_loo) - log10(EstMeanFlow)), na.rm=TRUE)

median(abs(log10(EstMeanFlow)     - log10(ObsMean)), na.rm=TRUE)
median(abs(log10(EstMeanFlow_loo) - log10(ObsMean)), na.rm=TRUE)


ObsMean[looL1model$site_no == "08167000"]*(1/0.3048^3)
EstMeanFlow[looL1model$site_no == "08167000"]*(1/0.3048^3)
EstMeanFlow_loo[looL1model$site_no == "08167000"]*(1/0.3048^3)

pdf("true_mean_comparison.pdf", useDingbats=FALSE, height=6, width=7.5)
par(mgp=c(3,0.5,0)) # going to tick to the inside, change some parameters
xlim <- add.log.axis(x=ObsMean); ylim <- add.log.axis(x=EstMeanFlow_loo)
plot(ObsMean, EstMeanFlow_loo, log="xy", pch=16, lwd=0.5, cex=1,
     xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i", xlim=xlim, ylim=ylim,
     col=rgb(1-as.numeric(looPPLOmodel$pplo == 0),0,as.numeric(looPPLOmodel$pplo == 0),.3))
abline(0,1)
add.log.axis(side=2,      tcl=0.8*abs(par()$tcl), two.sided=TRUE)
#add.log.axis(logs=c(1),   tcl=-0.5*abs(par()$tcl), side=2, two.sided=TRUE)
add.log.axis(logs=c(1),   tcl=+1.3*abs(par()$tcl), side=2, two.sided=TRUE)
add.log.axis(side=1,      tcl=0.8*abs(par()$tcl), two.sided=TRUE)
#add.log.axis(logs=c(1),   tcl=-0.5*abs(par()$tcl), side=1, two.sided=TRUE)
add.log.axis(logs=c(1),   tcl=+1.3*abs(par()$tcl), side=1, two.sided=TRUE)
add.log.axis(logs=c(1, 2, 3, 4, 6), side=2, make.labs=TRUE, las=1,
             label="GAM L1+PPLO LOO estimated decadal mean flow, in cms")
add.log.axis(logs=c(1, 2, 3, 4, 6), side=1, make.labs=TRUE, las=1,
             label="Observed decadal mean flow (zeros included), in cms")

points(ObsMean[        looL1model$site_no == "08167000"],
       EstMeanFlow[    looL1model$site_no == "08167000"], col=rgb(0,1,0,.8), pch=16, cex=1.25)
points(ObsMean[        looL1model$site_no == "08167000"],
       EstMeanFlow_loo[looL1model$site_no == "08167000"], col=rgb(0,1,0,.8), pch=1, lwd=2, cex=2)

legend(0.011,800, c("Equal value line",
                    "GAM L1 LOO streamflow for which no no-flow days were observed",
                    "GAM L1 LOO streamflow for which some no-flow days were observed",
                    "USGS streamgage 08167000 (GAM L1+PPLO estimated [non-LOO])",
                    "USGS streamgage 08167000 (GAM L1+PPLO LOO estimated)"),
          bty="n", lty=c(1,NA,NA,NA,NA), col=c(1,4,2,3,3), pch=c(NA,16,16,16,1),
          cex = 0.8, pt.cex=c(NA,1,1,1.25,2))
txt <- paste0("Abbreviations: GAM, generalized additive model;\n",
              "   L1, nonzero-flow mean; PPLO, decadal percentage of no flow;\n",
              "   LOO, leave-one-watershed-out, cross-validation of GAMs; and\n",
              "   cms, cubic meters per second")
text(.3, 0.03, txt, pos=4, cex=0.8)
par(mgp=c(3,1,0)) # restore defaults
dev.off()


sum(looPPLOmodel$pplo         >  0            ) # [1] 784
sum(looPPLOmodel$pplo         == 0            ) # [1] 2020
sum(looPPLOmodel$est_pplo     >  0            ) # [1] 540
sum(looPPLOmodel$est_pplo     == 0            ) # [1] 2264
sum(looPPLOmodel$loo_est_pplo >  0, na.rm=TRUE) # [1] 511
sum(looPPLOmodel$loo_est_pplo == 0, na.rm=TRUE) # [1] 2235


summary(looPPLOmodel$pplo        )
summary(looPPLOmodel$est_pplo    )
summary(looPPLOmodel$loo_est_pplo)



summary(looPPLOmodel$est_lwr_pplo > looPPLOmodel$flowtime |
        looPPLOmodel$est_upr_pplo < looPPLOmodel$flowtime )



n <- length(looL1model$L1)
o1 <- sum(looL1model$est_lwr_L1 > looL1model$L1 |
          looL1model$est_upr_L1 < looL1model$L1, na.rm=TRUE)
o2 <- sum(looL1model$loo_est_lwr_L1 > looL1model$L1 |
          looL1model$loo_est_upr_L1 < looL1model$L1, na.rm=TRUE)
cp1 <- round(1 - o1/n, 3)*100; cp2 <- round(1 - o2/n, 3)*100
message("95% coverage probability, L1: GAM=",cp1," and GAMloo=",cp2)

n <- length(looT2model$T2)
o1 <- sum(looT2model$est_lwr_T2 > looT2model$T2 |
          looT2model$est_upr_T2 < looT2model$T2, na.rm=TRUE)
o2 <- sum(looT2model$loo_est_lwr_T2 > looT2model$T2 |
          looT2model$loo_est_upr_T2 < looT2model$T2, na.rm=TRUE)
cp1 <- round(1 - o1/n, 3)*100; cp2 <- round(1 - o2/n, 3)*100
message("95% coverage probability, T2: GAM=",cp1," and GAMloo=",cp2)

n <- length(looT3model$T3)
o1 <- sum(looT3model$est_lwr_T3 > looT3model$T3 |
          looT3model$est_upr_T3 < looT3model$T3, na.rm=TRUE)
o2 <- sum(looT3model$loo_est_lwr_T3 > looT3model$T3 |
          looT3model$loo_est_upr_T3 < looT3model$T3, na.rm=TRUE)
cp1 <- round(1 - o1/n, 3)*100; cp2 <- round(1 - o2/n, 3)*100
message("95% coverage probability, T3: GAM=",cp1," and GAMloo=",cp2)

n <- length(looT4model$T4)
o1 <- sum(looT4model$est_lwr_T4 > looT4model$T4 |
          looT4model$est_upr_T4 < looT4model$T4, na.rm=TRUE)
o2 <- sum(looT4model$loo_est_lwr_T4 > looT4model$T4 |
          looT4model$loo_est_upr_T4 < looT4model$T4, na.rm=TRUE)
cp1 <- round(1 - o1/n, 3)*100; cp2 <- round(1 - o2/n, 3)*100
message("95% coverage probability, T4: GAM=",cp1," and GAMloo=",cp2)

n <- length(looT5model$T5)
o1 <- sum(looT5model$est_lwr_T5 > looT5model$T5 |
          looT5model$est_upr_T5 < looT5model$T5, na.rm=TRUE)
o2 <- sum(looT5model$loo_est_lwr_T5 > looT5model$T5 |
          looT5model$loo_est_upr_T5 < looT5model$T5, na.rm=TRUE)
cp1 <- round(1 - o1/n, 3)*100; cp2 <- round(1 - o2/n, 3)*100
message("95% coverage probability, T5: GAM=",cp1," and GAMloo=",cp2)
