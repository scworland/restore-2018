library(lmomco)
library(feather)

load("./Models.RData")
#looOverL1model <- read_feather("../../../results/gage/gam/all_gage_looest_overL1.feather")
looPPLOmodel <- read_feather("../../../results/gage/gam/all_gage_looest_pplo.feather")
looL1model   <- read_feather("../../../results/gage/gam/all_gage_looest_L1.feather"  )
looT2model   <- read_feather("../../../results/gage/gam/all_gage_looest_T2.feather"  )
looT3model   <- read_feather("../../../results/gage/gam/all_gage_looest_T3.feather"  )
looT4model   <- read_feather("../../../results/gage/gam/all_gage_looest_T4.feather"  )
looT5model   <- read_feather("../../../results/gage/gam/all_gage_looest_T5.feather"  )

chkBasins <- data.frame(site_no=looL1model$site_no, decade=looL1model$decade,
                        stringsAsFactors=FALSE)
areas <- sapply(1:length(chkBasins$site_no),
               function(i) { d <- chkBasins$decade[i]; s <- chkBasins$site_no[i]
                             DDo$basin_area[DDo$site_no == s & DDo$decade == d]  })
for(s in unique(DDo$site_no)) {
  for(d in DDo$decade[DDo$site_no == s]) {
    rt <- looL1model$L1[looL1model$site_no == s & looL1model$decade == d]
    if(length(rt) == 0) {
      message("Not found ",s," and ",d)
      # We don't want to find 02341500 for 5 decades as 02341505 has that record too
    } else {
      #message(s," and ",d," and ",length(rt))
    }
  }
}


EstMeanFlow <- looOverL1model$est_overL1
EstMeanFlow_loo <- looOverL1model$loo_est_overL1

ObsMean <- (1-looPPLOmodel$pplo)*looL1model$L1

median(abs(log10(EstMeanFlow_loo) - log10(EstMeanFlow)), na.rm=TRUE)
median(abs(log10(EstMeanFlow)     - log10(ObsMean)), na.rm=TRUE)
median(abs(log10(EstMeanFlow_loo) - log10(ObsMean)), na.rm=TRUE)

pdf("true_mean_comparison.pdf", useDingbats=FALSE, height=6, width=7.5)
par(mgp=c(3,0.5,0)) # going to tick to the inside, change some parameters
xlim <- add.log.axis(x=ObsMean); ylim <- add.log.axis(x=EstMeanFlow_loo)
xlim <- ylim <- c(6E-3, 1E3)
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

site <- "08167000"
#for(site in unique(DDo$site_no[DDo$ed_rch_zone == 1])) {
gcol <- 3
#if(site == "08155300" | site == "08155400" | site == "08156800") gcol <- 0
#site <- "08155300" # on the line, delist from Edwards?
#site <- "08155400" # on the line, delist from Edwards?
#site <- "08156800" # on the line, delist from Edwards?
#site <- "08181400" # close to the line
#site <- "08184000" # far from to the line
#site <- "08185000" # far from to the line
#site <- "08190500" # far from to the line
#site <- "08197500" # far from to the line
#site <- "08198500" # far from to the line
#site <- "08200700" # far from to the line
#site <- "08202700" # far from to the line
#print(mean(DDo$pplo[DDo$site_no == site]))

ObsMean[        looL1model$site_no == site]
EstMeanFlow[    looL1model$site_no == site]
EstMeanFlow_loo[looL1model$site_no == site]
ObsMean[        looL1model$site_no == site]*(1/0.3048^3)
EstMeanFlow[    looL1model$site_no == site]*(1/0.3048^3)
EstMeanFlow_loo[looL1model$site_no == site]*(1/0.3048^3)
green <- GISTools::add.alpha(rgb(0.5,1,0.5), 0.8)
points(ObsMean[        looL1model$site_no == site],
       EstMeanFlow[    looL1model$site_no == site], col=green, pch=16, cex=1.25)
points(ObsMean[        looL1model$site_no == site],
       EstMeanFlow_loo[looL1model$site_no == site], col=green, pch=1, lwd=2, cex=2)
#} # end for
legend(0.011,800, c("Equal value line",
                    "GAM L1 LOO streamflow for which no no-flow days were observed",
                    "GAM L1 LOO streamflow for which some no-flow days were observed",
                    paste0("USGS streamgage ",site," (GAM L1+PPLO estimated [non-LOO])"),
                    paste0("USGS streamgage ",site," (GAM L1+PPLO LOO estimated)")),
          bty="n", lty=c(1,NA,NA,NA,NA), col=c(1,4,2,rgb(0.5,1,0.5),rgb(0.5,1,0.5)), pch=c(NA,16,16,16,1),
          cex = 0.8, pt.cex=c(NA,1,1,1.25,2))
txt <- paste0("Abbreviations: GAM, generalized additive model;\n",
              "   L1, nonzero-flow mean; PPLO, decadal percentage of no flow;\n",
              "   LOO, leave-one-watershed-out, cross-validation of GAMs; and\n",
              "   cms, cubic meters per second")
text(.3, 0.03, txt, pos=4, cex=0.8)
par(mgp=c(3,1,0)) # restore defaults
dev.off()

# The coverage probability for the overall mean is too large
tf <- looOverL1model$in_model_L1 == 1 & looOverL1model$in_model_pplo
tmp <- looOverL1model[tf,]; tmpo <- ObsMean[tf]
n <- length(tmpo)
o1 <- sum(tmp$est_lwr_overL1     > tmpo | tmp$est_upr_overL1     < tmpo)
o2 <- sum(tmp$loo_est_lwr_overL1 > tmpo | tmp$loo_est_upr_overL1 < tmpo)
cp1 <- round(1 - o1/n, 3); cp2 <- round(1 - o2/n, 3)
message("95% coverage probability, overL1: GAM=",cp1," and GAMloo=",cp2)
corr <- 0.23
o1 <- sum((1+corr)*tmp$est_lwr_overL1     > tmpo | (1-corr)*tmp$est_upr_overL1     < tmpo)
o2 <- sum((1+corr)*tmp$loo_est_lwr_overL1 > tmpo | (1-corr)*tmp$loo_est_upr_overL1 < tmpo)
cp1 <- round(1 - o1/n, 3); cp2 <- round(1 - o2/n, 3)
message("Corrected 95% coverage probability, overL1: GAM=",cp1," and GAMloo=",cp2)

looOverL1model_corr <- looOverL1model
looOverL1model_corr$est_lwr_overL1     <- (1+corr)*looOverL1model_corr$est_lwr_overL1
looOverL1model_corr$est_upr_overL1     <- (1-corr)*looOverL1model_corr$est_upr_overL1
looOverL1model_corr$loo_est_lwr_overL1 <- (1+corr)*looOverL1model_corr$loo_est_lwr_overL1
looOverL1model_corr$loo_est_upr_overL1 <- (1-corr)*looOverL1model_corr$loo_est_upr_overL1
#write_feather(looOverL1model_corr, "all_gage_looest_overL1.feather")


tmp <- looL1model[looL1model$in_model_L1 == 1,]
n <- length(tmp$L1)
o1 <- sum(tmp$est_lwr_L1     > tmp$L1 | tmp$est_upr_L1     < tmp$L1)
o2 <- sum(tmp$loo_est_lwr_L1 > tmp$L1 | tmp$loo_est_upr_L1 < tmp$L1)
cp1 <- round(1 - o1/n, 3); cp2 <- round(1 - o2/n, 3)
message("95% coverage probability, L1: GAM=",cp1," and GAMloo=",cp2)

tmp <- looT2model[looT2model$in_model_T2 == 1,]
n <- length(tmp$T2)
o1 <- sum(tmp$est_lwr_T2     > tmp$T2 | tmp$est_upr_T2     < tmp$T2)
o2 <- sum(tmp$loo_est_lwr_T2 > tmp$T2 | tmp$loo_est_upr_T2 < tmp$T2)
cp1 <- round(1 - o1/n, 3); cp2 <- round(1 - o2/n, 3)
message("95% coverage probability, T2: GAM=",cp1," and GAMloo=",cp2)

tmp <- looT3model[looT3model$in_model_T3 == 1,]
n <- length(tmp$T3)
o1 <- sum(tmp$est_lwr_T3     > tmp$T3 | tmp$est_upr_T3     < tmp$T3)
o2 <- sum(tmp$loo_est_lwr_T3 > tmp$T3 | tmp$loo_est_upr_T3 < tmp$T3)
cp1 <- round(1 - o1/n, 3); cp2 <- round(1 - o2/n, 3)
message("95% coverage probability, T3: GAM=",cp1," and GAMloo=",cp2)

tmp <- looT4model[looT4model$in_model_T4 == 1,]
n <- length(tmp$T4)
o1 <- sum(tmp$est_lwr_T4     > tmp$T4 | tmp$est_upr_T4     < tmp$T4)
o2 <- sum(tmp$loo_est_lwr_T4 > tmp$T4 | tmp$loo_est_upr_T4 < tmp$T4)
cp1 <- round(1 - o1/n, 3); cp2 <- round(1 - o2/n, 3)
message("95% coverage probability, T4: GAM=",cp1," and GAMloo=",cp2)

tmp <- looT5model[looT5model$in_model_T5 == 1,]
n <- length(tmp$T5)
o1 <- sum(tmp$est_lwr_T5     > tmp$T5 | tmp$est_upr_T5     < tmp$T5)
o2 <- sum(tmp$loo_est_lwr_T5 > tmp$T5 | tmp$loo_est_upr_T5 < tmp$T5)
cp1 <- round(1 - o1/n, 3); cp2 <- round(1 - o2/n, 3)
message("95% coverage probability, T5: GAM=",cp1," and GAMloo=",cp2)





tmp <- looL1model[looL1model$in_model_L1 == 1,]
message("L1 RMSE: ",
round(sqrt(mean((log10(tmp$L1)-log10(tmp$est_L1))^2)), 3), " & ",
round(sqrt(mean((log10(tmp$L1)-log10(tmp$loo_est_L1))^2)), 3)
)
tmp <- looT2model[looT2model$in_model_T2 == 1,]
message("T2 RMSE: ",
round(sqrt(mean((tmp$T2-tmp$est_T2)^2, na.rm=TRUE)), 3), " & ",
round(sqrt(mean((tmp$T2-tmp$loo_est_T2)^2, na.rm=TRUE)), 3)
)
tmp <- looT3model[looT3model$in_model_T3 == 1,]
message("T3 RMSE: ",
round(sqrt(mean((tmp$T3-tmp$est_T3)^2, na.rm=TRUE)), 3), " & ",
round(sqrt(mean((tmp$T3-tmp$loo_est_T3)^2, na.rm=TRUE)), 3)
)
tmp <- looT4model[looT4model$in_model_T4 == 1,]
message("T4 RMSE: ",
round(sqrt(mean((tmp$T4-tmp$est_T4)^2, na.rm=TRUE)), 3), " & ",
round(sqrt(mean((tmp$T4-tmp$loo_est_T4)^2, na.rm=TRUE)), 3)
)
tmp <- looT5model[looT5model$in_model_T5 == 1,]
message("T5 RMSE: ",
round(sqrt(mean((tmp$T5-tmp$est_T5)^2, na.rm=TRUE)), 3), " & ",
round(sqrt(mean((tmp$T5-tmp$loo_est_T5)^2, na.rm=TRUE)), 3)
)


tmp <- looPPLOmodel[looPPLOmodel$in_model_pplo == 1,]
n <- length(tmp$pplo); print(n)
sum(tmp$pplo > 0)/n
sum(tmp$est_pplo > 0)/n
sum(tmp$loo_est_pplo > 0)/n
(sum(tmp$pplo > 0 & tmp$est_pplo > 0) + sum(tmp$pplo == 0 & tmp$est_pplo == 0))/n
(sum(tmp$pplo > 0 & tmp$loo_est_pplo > 0) + sum(tmp$pplo == 0 & tmp$loo_est_pplo == 0))/n

mean(abs(tmp$est_pplo - tmp$pplo) <= 0   )
mean(abs(tmp$est_pplo - tmp$pplo) <= 0.02)
mean(abs(tmp$est_pplo - tmp$pplo) <= 0.05)
mean(abs(tmp$est_pplo - tmp$pplo) <= 0.10)

mean(abs(tmp$loo_est_pplo - tmp$pplo) <= 0   )
mean(abs(tmp$loo_est_pplo - tmp$pplo) <= 0.02)
mean(abs(tmp$loo_est_pplo - tmp$pplo) <= 0.05)
mean(abs(tmp$loo_est_pplo - tmp$pplo) <= 0.10)


message("PPLO RMSE (flowtime): ",
round(sqrt(mean((tmp$flowtime-tmp$est_flowtime)^2)), 3), " & ",
round(sqrt(mean((tmp$flowtime-tmp$loo_est_flowtime)^2)), 3)
)


sum(tmp$pplo         >  0) # [1] 748
sum(tmp$pplo         == 0) # [1] 2002
sum(tmp$est_pplo     >  0) # [1] 510
sum(tmp$est_pplo     == 0) # [1] 2240
sum(tmp$loo_est_pplo >  0) # [1] 518
sum(tmp$loo_est_pplo == 0) # [1] 2232


summary(tmp$pplo        )
summary(tmp$est_pplo    )
summary(tmp$loo_est_pplo)


n <- length(tmp$flowtime)
o1 <- sum(tmp$est_lwr_flowtime > tmp$flowtime |
          tmp$est_upr_flowtime < tmp$flowtime)
o2 <- sum(tmp$loo_est_lwr_flowtime > tmp$flowtime |
          tmp$loo_est_upr_flowtime < tmp$flowtime)
cp1 <- round(1 - o1/n, 3); cp2 <- round(1 - o2/n, 3)
message("95% coverage probability, PPLO: GAM=",cp1," and GAMloo=",cp2)



message("PPLO RMSE (pplo): ",
round(sqrt(mean((tmp$pplo - tmp$est_pplo)^2)), digits=3), " & ",
round(sqrt(mean((tmp$pplo - tmp$loo_est_pplo)^2)), digits=3))

jtmpa <- tmp[tmp$est_pplo     != 0,] # Kroll and Stedinger (1999)
jtmpb <- tmp[tmp$loo_est_pplo != 0,] # Kroll and Stedinger (1999)
message("PPLO RMSE (pplo [KS1999]): ",
round(sqrt(mean((jtmpa$pplo - jtmpa$est_pplo)^2)), digits=3), " & ",
round(sqrt(mean((jtmpb$pplo - jtmpb$loo_est_pplo)^2)) , digits=3))
rm(jtmpa, jtmpb)


whole.correct <- (sum(tmp$pplo  > 0 & tmp$est_pplo  > 0) +
                  sum(tmp$pplo == 0 & tmp$est_pplo == 0))/n
loo.correct   <- (sum(tmp$pplo  > 0 & tmp$loo_est_pplo  > 0) +
                  sum(tmp$pplo == 0 & tmp$loo_est_pplo == 0))/n

message("PPLO -- lengths(site_no, comid, huc12): ",
length(unique(looPPLOmodel$site_no)), ", ",
length(unique(looPPLOmodel$comid)), ", ",
length(unique(looPPLOmodel$huc12)))

message("L1 -- lengths(site_no, comid, huc12): ",
length(unique(looL1model$site_no)), ", ",
length(unique(looL1model$comid)), ", ",
length(unique(looL1model$huc12)))

message("T2 -- lengths(site_no, comid, huc12): ",
length(unique(looT2model$site_no)), ", ",
length(unique(looT2model$comid)), ", ",
length(unique(looT2model$huc12)))

message("T3 -- lengths(site_no, comid, huc12): ",
length(unique(looT3model$site_no)), ", ",
length(unique(looT3model$comid)), ", ",
length(unique(looT3model$huc12)))

message("T4 -- lengths(site_no, comid, huc12): ",
length(unique(looT4model$site_no)), ", ",
length(unique(looT4model$comid)), ", ",
length(unique(looT4model$huc12)))

message("T5 -- lengths(site_no, comid, huc12): ",
length(unique(looT5model$site_no)), ", ",
length(unique(looT5model$comid)), ", ",
length(unique(looT5model$huc12)))

