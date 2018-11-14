library(lmomco)
library(feather)

hPPLO <- read_feather("../../../results/huc12/gam/all_huc12_pplo.feather")
hL1   <- read_feather("../../../results/huc12/gam/all_huc12_L1.feather")
hT2   <- read_feather("../../../results/huc12/gam/all_huc12_T2.feather")
hT3   <- read_feather("../../../results/huc12/gam/all_huc12_T3.feather")
hT4   <- read_feather("../../../results/huc12/gam/all_huc12_T4.feather")
hT5   <- read_feather("../../../results/huc12/gam/all_huc12_T5.feather")

G   <- read_feather("../../../data/gage/all_gage_data.feather")

fdc <- as.data.frame(G[,c(11:37)])
ff <- as.numeric(gsub("f","",names(fdc)))/100

tmp <- hT2; tmp$est_lwr_T2 <- NULL; tmp$est_upr_T2 <- NULL
tmp$est_T2 <- NULL; tmp$rse_T2 <- NULL; tmp$se.fit_T2 <- NULL
tmp$delta_est_T2 <- NULL
jmp <- as.list(G[1,c(11:37)])
jmp$f0.02 <- rep(NA, length(tmp$comid))
jmp$f0.05 <- rep(NA, length(tmp$comid))
jmp$f0.1 <- rep(NA, length(tmp$comid))
jmp$f0.2 <- rep(NA, length(tmp$comid))
jmp$f0.5 <- rep(NA, length(tmp$comid))
jmp$f01 <- rep(NA, length(tmp$comid))
jmp$f02 <- rep(NA, length(tmp$comid))
jmp$f05 <- rep(NA, length(tmp$comid))
jmp$f10 <- rep(NA, length(tmp$comid))
jmp$f20 <- rep(NA, length(tmp$comid))
jmp$f25 <- rep(NA, length(tmp$comid))
jmp$f30 <- rep(NA, length(tmp$comid))
jmp$f40 <- rep(NA, length(tmp$comid))
jmp$f50 <- rep(NA, length(tmp$comid))
jmp$f60 <- rep(NA, length(tmp$comid))
jmp$f70 <- rep(NA, length(tmp$comid))
jmp$f75 <- rep(NA, length(tmp$comid))
jmp$f80 <- rep(NA, length(tmp$comid))
jmp$f90 <- rep(NA, length(tmp$comid))
jmp$f95 <- rep(NA, length(tmp$comid))
jmp$f98 <- rep(NA, length(tmp$comid))
jmp$f99 <- rep(NA, length(tmp$comid))
jmp$f99.5 <- rep(NA, length(tmp$comid))
jmp$f99.8 <- rep(NA, length(tmp$comid))
jmp$f99.9 <- rep(NA, length(tmp$comid))
jmp$f99.95 <- rep(NA, length(tmp$comid))
jmp$f99.98 <- rep(NA, length(tmp$comid))
jmp <- as.data.frame(jmp)
tmp <- cbind(tmp,jmp)

PHgno <- PHkap <- PHaep4 <- tmp
PHkap$snapped <- PHaep4$snapped <- FALSE
PHgno$snapped <- "not applicable"
PHgno$est_pplo <- PHkap$est_pplo <- PHaep4$est_pplo <- NA
PHgno$overall_mean <- PHkap$overall_mean <- PHaep4$overall_mean <- NA

mc <- 1E7; ru <- runif(mc)
n <- length(PHaep4$comid)
ix <- 1:n
#ix[PHaep4$comid == "10622893" & PHaep4$decade == "2000"]
#ix <- ix[is.na(PHaep4$overall_mean)]

for(i in ix) {
  message("Working on ",PHaep4$comid[i]," and ",PHaep4$decade[i],
           ", ",i,"(",n,")")
  suppressWarnings(lmr <- vec2lmom(c(hL1$est_L1[i]*hL1$bias_corr[i],
   hT2$est_T2[i], hT3$est_T3[i], hT4$est_T4[i], hT5$est_T5[i]), lscale=FALSE))
  if(! are.lmom.valid(lmr)) {
    #message("   invalid L-moments for ",PHaep4$comid[i]," and ",PHaep4$decade[i])
    if(lmr$ratios[4] <  (5 * lmr$ratios[3]^2 - 1)/4) {
       lmr$ratios[4] <- (5 * lmr$ratios[3]^2 - 1)/4
    }
    if(lmr$ratios[3] < -1) lmr$ratios[3] <- -0.99
    if(lmr$ratios[3] > +1) lmr$ratios[3] <- +0.99
    if(lmr$ratios[4] < -1) lmr$ratios[4] <- -0.99
    if(lmr$ratios[4] > +1) lmr$ratios[4] <- +0.99
    if(lmr$ratios[5] < -1) lmr$ratios[5] <- -0.99
    if(lmr$ratios[5] > +1) lmr$ratios[5] <- +0.99
  }
  my.pplo <- hPPLO$est_pplo[i]
  PHgno$est_pplo[i] <- PHkap$est_pplo[i] <- PHaep4$est_pplo[i] <- my.pplo

  aep4 <- paraep4(lmr, snap.tau4=TRUE)
  if(length(aep4$message) != 0) {
    message(" below AEP4 snapped to lower Tau4 bounds")
    PHaep4$snapped[i] <- TRUE
  }
  gno  <- pargno(lmr)
  kap <- parkap(lmr, snap.tau4=FALSE)
  if(kap$ifail > 0) {
     message(" above KAP snapped to upper Tau4 bounds")
     PHkap$snapped[i] <- TRUE
  }
  kap <- parkap(lmr, snap.tau4=TRUE)

  if(! is.null(aep4)) {
     qua1 <- qlmomco(f2flo(ff, pp=my.pplo), aep4)
     qua1[qua1 < 0] <- 0
     qua1 <- c(rep(0,27-length(qua1)),qua1)
     for(j in 6:32) PHaep4[i,j] <- qua1[j-5]
  } else {
     qua1 <- rep(NA,27)
     message("aep4 error")
  }

  if(! is.null(gno)) {
     qua2 <- qlmomco(f2flo(ff, pp=my.pplo), gno)
     qua2[qua2 < 0] <- 0
     qua2 <- c(rep(0,27-length(qua2)),qua2)
     for(j in 6:32) PHgno[i,j] <- qua2[j-5]
  } else {
     qua2 <- rep(NA,27)
     message("gno error")
  }
  if(kap$ifail == 0) {
    qua3 <- qlmomco(f2flo(ff, pp=my.pplo), kap)
    qua3[qua3 < 0] <- 0
    qua3 <- c(rep(0,27-length(qua3)),qua3)
    for(j in 6:32) PHkap[i,j] <- qua3[j-5]
  } else {
    qua3 <- rep(NA,27)
    message("kappa error")
  }
  qua1[is.na(qua1)] <- 0; qua2[is.na(qua2)] <- 0; qua3[is.na(qua3)] <- 0
  Mu1 <- function(f) approx(ff,qua1,xout=f,rule=2)$y
  mu1 <- NULL
  try(mu1 <- integrate(Mu1, lower=0, upper=1)$value, silent=TRUE)
  if(is.null(mu1)) {
     mu1 <- sum(approx(ff,qua1,xout=ru,rule=2)$y)/mc
  }
  PHaep4$overall_mean[i] <- mu1

  Mu2 <- function(f) approx(ff,qua2,xout=f,rule=2)$y
  mu2 <- NULL
  try(mu2 <- integrate(Mu2, lower=0, upper=1)$value, silent=TRUE)
  if(is.null(mu2)) {
     mu2 <- sum(approx(ff,qua2,xout=ru,rule=2)$y)/mc
  }
  PHgno$overall_mean[i] <- mu2

  Mu3 <- function(f) approx(ff,qua3,xout=f,rule=2)$y
  mu3 <- NULL
  try(mu3 <- integrate(Mu3, lower=0, upper=1)$value, silent=TRUE)
  if(is.null(mu3)) {
     mu3 <- sum(approx(ff,qua3,xout=ru,rule=2)$y)/mc
  }
  PHkap$overall_mean[i] <- mu3

  message("   *** mus ", mu1, "(aep4), ", mu2, "(gno), ", mu3, "(kap)")
}

PHaep4$overall_mean[PHaep4$overall_mean == 0] <- NA
PHgno$overall_mean[PHgno$overall_mean == 0] <- NA
PHkap$overall_mean[PHkap$overall_mean == 0] <- NA

for(i in 1:length(PGaep4$comid)) {
  for(j in 7:33) {
    if(PGaep4[i,j] == -9999) PGaep4[i,j] <- NA
    if(PGgno[i,j]  == -9999) PGgno[i,j] <- NA
    if(PGkap[i,j]  == -9999) PGkap[i,j] <- NA
  }
}



write_feather(PHaep4, "all_huc_fdclmr_aep.feather")
write_feather(PHgno,  "all_huc_fdclmr_gno.feather")
write_feather(PHkap,  "all_huc_fdclmr_kap.feather")

