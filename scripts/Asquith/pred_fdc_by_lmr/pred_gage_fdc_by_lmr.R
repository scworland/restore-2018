library(lmomco)
library(feather)
gPPLO <- read_feather("../../../results/gage/gam/all_gage_looest_pplo.feather")
gL1   <- read_feather("../../../results/gage/gam/all_gage_looest_L1.feather")
gT2   <- read_feather("../../../results/gage/gam/all_gage_looest_T2.feather")
gT3   <- read_feather("../../../results/gage/gam/all_gage_looest_T3.feather")
gT4   <- read_feather("../../../results/gage/gam/all_gage_looest_T4.feather")
gT5   <- read_feather("../../../results/gage/gam/all_gage_looest_T5.feather")

G   <- read_feather("../../../data/gage/all_gage_data.feather")

fdc <- as.data.frame(G[,c(11:37)])
ff <- as.numeric(gsub("f","",names(fdc)))/100


PGaep4 <- G[,c(1:6,11:37)]
PGaep4[,7:33] <- -9999
PGgno <- PGkap <- PGaep4
PGkap$snapped <- PGaep4$snapped <- FALSE
PGgno$snapped <- "not applicable"
PGgno$est_pplo <- PGkap$est_pplo <- PGaep4$est_pplo <- NA
PGgno$overall_mean <- PGkap$overall_mean <- PGaep4$overall_mean <- NA

mc <- 1E7; ru <- runif(mc)
n <- length(PGaep4$comid)
ix <- 1:n


for(i in ix) {
  if(length(gL1$site_no[gL1$site_no == PGaep4$site_no[i]]) == 0) {
    message("**** Not found site=",PGaep4$site_no[i]) # 02341500
    # same records as 02341505 but later has one more decade so
    # 02341500 is NOT represented in predictions of the L-moments and thus
    # not the FDC.
  }
  message("Working on ",PGaep4$site_no[i]," and ",PGaep4$decade[i],
           ", ",i,"(",n,")")
  if(is.na(gL1$est_L1[i])) next # 08211000 and 1980, 2800(2804)
  suppressWarnings(lmr <- vec2lmom(c(gL1$est_L1[i]*gL1$bias_corr[i],
   gT2$est_T2[i], gT3$est_T3[i], gT4$est_T4[i], gT5$est_T5[i]), lscale=FALSE))
  if(! are.lmom.valid(lmr)) {
    #message("   invalid L-moments for ",PGaep4$site_no[i]," and ",PGaep4$decade[i])
    if(lmr$ratios[4] < (5 * lmr$ratios[3]^2 - 1)/4) {
       lmr$ratios[4] <- (5 * lmr$ratios[3]^2 - 1)/4
    }
    if(lmr$ratios[3] < -1) lmr$ratios[3] <- -0.99
    if(lmr$ratios[3] > +1) lmr$ratios[3] <- +0.99
    if(lmr$ratios[4] < -1) lmr$ratios[4] <- -0.99
    if(lmr$ratios[4] > +1) lmr$ratios[4] <- +0.99
    if(lmr$ratios[5] < -1) lmr$ratios[5] <- -0.99
    if(lmr$ratios[5] > +1) lmr$ratios[5] <- +0.99
  }
  my.pplo <- gPPLO$est_pplo[i]
  PGgno$est_pplo[i] <- PGkap$est_pplo[i] <- PGaep4$est_pplo[i] <- my.pplo

  aep4 <- paraep4(lmr, snap.tau4=TRUE)
  if(length(aep4$message) != 0) {
    #message(" below AEP4 snapped to lower Tau4 bounds")
    PGaep4$snapped[i] <- TRUE
  }
  gno  <- pargno(lmr)
  kap <- parkap(lmr, snap.tau4=FALSE)
  if(kap$ifail > 0) {
     #message(" above KAP snapped to upper Tau4 bounds")
     PGkap$snapped[i] <- TRUE
  }
  kap <- parkap(lmr, snap.tau4=TRUE)

  if(! is.null(aep4)) {
     qua1 <- qlmomco(f2flo(ff, pp=my.pplo), aep4)
     qua1[qua1 < 0] <- 0
     qua1 <- c(rep(0,27-length(qua1)),qua1)
     for(j in 7:33) PGaep4[i,j] <- qua1[j-6]
  } else {
     qua1 <- rep(NA,27)
     #message("aep4 error")
  }

  if(! is.null(gno)) {
     qua2 <- qlmomco(f2flo(ff, pp=my.pplo), gno)
     qua2[qua2 < 0] <- 0
     qua2 <- c(rep(0,27-length(qua2)),qua2)
     for(j in 7:33) PGgno[i,j] <- qua2[j-6]
  } else {
     qua2 <- rep(NA,27)
     #message("gno error")
  }

  if(kap$ifail == 0) {
    qua3 <- qlmomco(f2flo(ff, pp=my.pplo), kap)
    qua3[qua3 < 0] <- 0
    qua3 <- c(rep(0,27-length(qua3)),qua3)
    for(j in 7:33) PGkap[i,j] <- qua3[j-6]
  } else {
    qua3 <- rep(NA,27)
    #message("kappa error")
  }

  qua1[is.na(qua1)] <- 0; qua2[is.na(qua2)] <- 0; qua3[is.na(qua3)] <- 0
  Mu1 <- function(f) approx(ff,qua1,xout=f,rule=2)$y
  mu1 <- NULL
  try(mu1 <- integrate(Mu1, lower=0, upper=1)$value, silent=TRUE)
  if(is.null(mu1)) {
     mu1 <- sum(approx(ff,qua1,xout=ru,rule=2)$y)/mc
  }
  PGaep4$overall_mean[i] <- mu1

  Mu2 <- function(f) approx(ff,qua2,xout=f,rule=2)$y
  mu2 <- NULL
  try(mu2 <- integrate(Mu2, lower=0, upper=1)$value, silent=TRUE)
  if(is.null(mu2)) {
     mu2 <- sum(approx(ff,qua2,xout=ru,rule=2)$y)/mc
  }
  PGgno$overall_mean[i] <- mu2

  Mu3 <- function(f) approx(ff,qua3,xout=f,rule=2)$y
  mu3 <- NULL
  try(mu3 <- integrate(Mu3, lower=0, upper=1)$value, silent=TRUE)
  if(is.null(mu3)) {
     mu3 <- sum(approx(ff,qua3,xout=ru,rule=2)$y)/mc
  }
  PGkap$overall_mean[i] <- mu3

  message("   *** mus ", mu1, "(aep4), ", mu2, "(gno), ", mu3, "(kap)")
}

PGaep4$snapped <- as.numeric(PGaep4$snapped)
PGgno$snapped <- as.numeric(PGaep4$snapped)
PGkap$snapped <- as.numeric(PGaep4$snapped)

PGaep4$overall_mean[PGaep4$overall_mean == 0] <- NA
PGgno$overall_mean[PGgno$overall_mean == 0] <- NA
PGkap$overall_mean[PGkap$overall_mean == 0] <- NA


for(i in 1:length(PGaep4$comid)) {
  for(j in 7:33) {
    if(PGaep4[i,j] == -9999) PGaep4[i,j] <- NA
    if(PGgno[i,j]  == -9999) PGgno[i,j] <- NA
    if(PGkap[i,j]  == -9999) PGkap[i,j] <- NA
  }
}


write_feather(PGaep4, "all_gage_fdclmr_aep.feather")
write_feather(PGgno,  "all_gage_fdclmr_gno.feather")
write_feather(PGkap,  "all_gage_fdclmr_kap.feather")

