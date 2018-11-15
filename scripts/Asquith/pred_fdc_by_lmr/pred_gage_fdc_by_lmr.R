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
PGgno$overall_mean_by_pploubL1 <- PGkap$overall_mean_by_pploubL1 <- PGaep4$overall_mean_by_pploubL1 <- NA
PGgno$overall_mean_by_dist <- PGkap$overall_mean_by_dist <- PGaep4$overall_mean_by_dist <- NA
mc <- 1E7; ru <- runif(mc)
n <- length(PGaep4$comid)
ix <- 1:n

#plot(1,1, type="n", xlim=c(.01,1000), ylim=c(0.01,1000), log="xy")
for(i in ix) {
  site <- PGaep4$site_no[i]; decade <- PGaep4$decade[i]
  if(length(gL1$site_no[gL1$site_no == site]) == 0) {
    message("**** Not found site=",site) # 02341500
    # same records as 02341505 but later has one more decade so
    # 02341500 is NOT represented in predictions of the L-moments and thus
    # not the FDC.
    next
  }
  message("Working on ",site," and ",decade,", ",i,"(",n,")")

  # 08211000 and 1980, 2800(2804)
  if(is.na(gL1$est_L1[gL1$site_no == site & gL1$decade == decade])) next

  my.unbias <- gL1$est_L1[   gL1$site_no == site & gL1$decade == decade]*
               gL1$bias_corr[gL1$site_no == site & gL1$decade == decade]

  suppressWarnings(
     lmr <- vec2lmom(
     c(my.unbias,
       gT2$est_T2[   gT2$site_no == site & gT2$decade == decade],
       gT3$est_T3[   gT3$site_no == site & gT3$decade == decade],
       gT4$est_T4[   gT4$site_no == site & gT4$decade == decade],
       gT5$est_T5[   gT5$site_no == site & gT5$decade == decade]), lscale=FALSE))
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
  my.pplo <- gPPLO$est_pplo[gPPLO$site_no == site & gPPLO$decade == decade]
  PGgno$est_pplo[i] <- PGkap$est_pplo[i] <- PGaep4$est_pplo[i] <- my.pplo
  PGgno$overall_mean_by_pploubL1[i]  <- (1-my.pplo)*my.unbias
  PGkap$overall_mean_by_pploubL1[i]  <- (1-my.pplo)*my.unbias
  PGaep4$overall_mean_by_pploubL1[i] <- (1-my.pplo)*my.unbias

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
  PGaep4$overall_mean_by_dist[i] <- mu1

  Mu2 <- function(f) approx(ff,qua2,xout=f,rule=2)$y
  mu2 <- NULL
  try(mu2 <- integrate(Mu2, lower=0, upper=1)$value, silent=TRUE)
  if(is.null(mu2)) {
     mu2 <- sum(approx(ff,qua2,xout=ru,rule=2)$y)/mc
  }
  PGgno$overall_mean_by_dist[i] <- mu2

  Mu3 <- function(f) approx(ff,qua3,xout=f,rule=2)$y
  mu3 <- NULL
  try(mu3 <- integrate(Mu3, lower=0, upper=1)$value, silent=TRUE)
  if(is.null(mu3)) {
     mu3 <- sum(approx(ff,qua3,xout=ru,rule=2)$y)/mc
  }
  PGkap$overall_mean_by_dist[i] <- mu3

  mu <- (1-gPPLO$est_pplo[i])*gL1$est_L1[i]
  #message("   *** mus ", round(mu, digits=4), "(lmr), ", round(mu1, digits=4),
  #        "(aep4), ",   round(mu2, digits=4), "(gno), ", round(mu3, digits=4), "(kap)")
  #points(mu, mu1, lwd=0.5, col=2, cex=0.50)
  #points(mu, mu2, lwd=0.5, col=3, cex=1.00)
  #points(mu, mu3, lwd=0.5, col=4, cex=1.50)
}
est_lmr_overall <- PGaep4$overall_mean_by_pploubL1

PGaep4[is.na(PGaep4$site_no),]
PGgno[is.na(PGgno$site_no),]
PGkap[is.na(PGkap$site_no),]

PGaep4$snapped <- as.numeric(PGaep4$snapped)
PGgno$snapped <- as.numeric(PGaep4$snapped)
PGkap$snapped <- as.numeric(PGaep4$snapped)

PGaep4$overall_mean_by_dist[! is.na(PGaep4$overall_mean_by_dist) & PGaep4$overall_mean_by_dist == 0] <- NA
PGgno$overall_mean_by_dist[! is.na(PGgno$overall_mean_by_dist) & PGgno$overall_mean_by_dist == 0] <- NA
PGkap$overall_mean_by_dist[! is.na(PGkap$overall_mean_by_dist) & PGkap$overall_mean_by_dist == 0] <- NA

plot(overall_mean_by_dist, PGaep4$overall_mean_by_dist, log="xy")
plot(overall_mean_by_dist, PGgno$overall_mean_by_dist, log="xy")
plot(overall_mean_by_dist, PGkap$overall_mean_by_dist, log="xy")

del1 <- log10(PGaep4$overall_mean_by_dist)-log10(est_lmr_overall)
del2 <- log10(PGgno$overall_mean_by_dist)-log10(est_lmr_overall)
del3 <- log10(PGkap$overall_mean_by_dist)-log10(est_lmr_overall)
boxplot(del1)
boxplot(del2)
boxplot(del3)

SA <- PGaep4[! is.na(PGaep4$overall_mean_by_dist) & abs(del1) > 0.1,]
SG <- PGgno[ ! is.na(PGgno$overall_mean_by_dist ) & abs(del2) > 0.1,]
SK <- PGkap[ ! is.na(PGkap$overall_mean_by_dist ) & abs(del3) > 0.1,]

#plot(ff,SK[2,7:33], type="l")
#lines(ff,G[G$site_no == "08194200" & G$decade == 1990,11:37], col=2)

for(s in unique(SA$site_no)) {
  tmp <- SA[SA$site_no == s,]
  for(d in tmp$decade) {
     message("AEP4: resetting quantiles and overall mean to NA for ",s," and ",d,
             " because means don't align within 0.1 log10")
     for(i in 7:33) {
        PGaep4[PGaep4$site_no == s & PGaep4$decade == d, i] <- NA
     }
     PGaep4$overall_mean_by_dist[PGaep4$site_no == s & PGaep4$decade == d] <- NA
  }
}
for(s in unique(SG$site_no)) {
  tmp <- SG[SG$site_no == s,]
  for(d in tmp$decade) {
     message("GNO: resetting quantiles and overall mean to NA for ",s," and ",d,
             " because means don't align within 0.1 log10")
     for(i in 7:33) {
        PGgno[PGgno$site_no == s & PGgno$decade == d, i] <- NA
     }
     PGgno$overall_mean_by_dist[PGgno$site_no == s & PGgno$decade == d] <- NA
  }
}
for(s in unique(SK$site_no)) {
  tmp <- SK[SK$site_no == s,]
  for(d in tmp$decade) {
     message("KAP: resetting quantiles and overall mean to NA for ",s," and ",d,
             " because means don't align within 0.1 log10")
     for(i in 7:33) {
        PGkap[PGkap$site_no == s & PGkap$decade == d, i] <- NA
     }
     PGkap$overall_mean_by_dist[PGkap$site_no == s & PGkap$decade == d] <- NA
  }
}

plot(est_lmr_overall, PGaep4$overall_mean_by_dist, log="xy")
plot(est_lmr_overall, PGgno$overall_mean_by_dist,  log="xy")
plot(est_lmr_overall, PGkap$overall_mean_by_dist,  log="xy")

PGaep4[is.na(PGaep4$site_no),]
PGgno[ is.na(PGgno$site_no), ]
PGkap[ is.na(PGkap$site_no), ]


# this is a type of leak detection, we don't want any -9999 left over
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

PGaep4 <- read_feather("all_gage_fdclmr_aep.feather")
PGgno  <- read_feather("all_gage_fdclmr_gno.feather")
PGkap  <- read_feather("all_gage_fdclmr_kap.feather")

PGaep4 <- PGaep4[PGaep4$site_no != "02341500",]
PGgno  <- PGgno [PGgno$site_no  != "02341500",]
PGkap  <- PGkap[ PGkap$site_no  != "02341500",]
Gtmp <- G[G$site_no != "02341500",]

lmr_overall <- (1-Gtmp$pplo)*Gtmp$L1
plot(lmr_overall, est_lmr_overall, log="xy")


