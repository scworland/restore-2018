library(lmomco)
library(feather)
library(sp)
LATLONG <- paste0("+proj=longlat +ellps=GRS80 ",
                  "+datum=NAD83 +no_defs +towgs84=0,0,0")
LATLONG <- sp::CRS(LATLONG)
ALBEA <- paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 ",
                "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
ALBEA <- sp::CRS(ALBEA)

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
PHgno$est_L1 <- PHkap$est_L1 <- PHaep4$est_L1 <- NA
PHgno$overall_mean_by_pploubL1 <- PHkap$overall_mean_by_pploubL1 <- PHaep4$overall_mean_by_pploubL1 <- NA
PHgno$overall_mean_by_dist <- PHkap$overall_mean_by_dist <- PHaep4$overall_mean_by_dist <- NA

mc <- 1E7; ru <- runif(mc)
n <- length(PHaep4$comid)
ix <- 1:n
#ix[PHaep4$comid == "10622893" & PHaep4$decade == "2000"]
#ix <- ix[is.na(PHaep4$overall_mean)]

for(i in ix) {
  site <- PHaep4$site_no[i]; decade <- PHaep4$decade[i]
  message("Working on ",site," and ",decade, ", ",i,"(",n,")")
  my.unbias <- hL1$est_L1[i]*hL1$bias_corr[i]
  suppressWarnings(lmr <- vec2lmom(c(my.unbias,
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
  PHgno$overall_mean_by_pploubL1[i]  <- (1-my.pplo)*my.unbias
  PHkap$overall_mean_by_pploubL1[i]  <- (1-my.pplo)*my.unbias
  PHaep4$overall_mean_by_pploubL1[i] <- (1-my.pplo)*my.unbias

  aep4 <- paraep4(lmr, snap.tau4=TRUE)
  if(length(aep4$message) != 0) {
    #message(" below AEP4 snapped to lower Tau4 bounds")
    PHaep4$snapped[i] <- TRUE
  }
  gno  <- pargno(lmr)
  kap <- parkap(lmr, snap.tau4=FALSE)
  if(kap$ifail > 0) {
     #message(" above KAP snapped to upper Tau4 bounds")
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
     #message("aep4 error")
  }

  if(! is.null(gno)) {
     qua2 <- qlmomco(f2flo(ff, pp=my.pplo), gno)
     qua2[qua2 < 0] <- 0
     qua2 <- c(rep(0,27-length(qua2)),qua2)
     for(j in 6:32) PHgno[i,j] <- qua2[j-5]
  } else {
     qua2 <- rep(NA,27)
     #message("gno error")
  }
  if(kap$ifail == 0) {
    qua3 <- qlmomco(f2flo(ff, pp=my.pplo), kap)
    qua3[qua3 < 0] <- 0
    qua3 <- c(rep(0,27-length(qua3)),qua3)
    for(j in 6:32) PHkap[i,j] <- qua3[j-5]
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
  PHaep4$overall_mean_by_dist[i] <- mu1

  Mu2 <- function(f) approx(ff,qua2,xout=f,rule=2)$y
  mu2 <- NULL
  try(mu2 <- integrate(Mu2, lower=0, upper=1)$value, silent=TRUE)
  if(is.null(mu2)) {
     mu2 <- sum(approx(ff,qua2,xout=ru,rule=2)$y)/mc
  }
  PHgno$overall_mean_by_dist[i] <- mu2

  Mu3 <- function(f) approx(ff,qua3,xout=f,rule=2)$y
  mu3 <- NULL
  try(mu3 <- integrate(Mu3, lower=0, upper=1)$value, silent=TRUE)
  if(is.null(mu3)) {
     mu3 <- sum(approx(ff,qua3,xout=ru,rule=2)$y)/mc
  }
  PHkap$overall_mean_by_dist[i] <- mu3

  mu <- (1-hPPLO$est_pplo[i])*hL1$est_L1[i]
  #message("   *** mus ", round(mu1, digits=4), "(aep4), ", round(mu2, digits=4), "(gno), ", round(mu3, digits=4), "(kap)")
}
est_lmr_overall <- PHaep4$overall_mean_by_pploubL1

write_feather(PHaep4, "ORGall_huc_fdclmr_aep.feather")
write_feather(PHgno,  "ORGall_huc_fdclmr_gno.feather")
write_feather(PHkap,  "ORGall_huc_fdclmr_kap.feather")

PHaep4 <- read_feather("ORGall_huc_fdclmr_aep.feather")
PHgno  <- read_feather("ORGall_huc_fdclmr_gno.feather")
PHkap  <- read_feather("ORGall_huc_fdclmr_kap.feather")


PHaep4[is.na(PHaep4$comid),]
PHgno[ is.na(PHaep4$comid),]
PHkap[ is.na(PHaep4$comid),]

PHaep4$snapped <- as.numeric(PHaep4$snapped)
#PHgno$snapped <- as.numeric(PHgno$snapped) # character "not applicable"
PHkap$snapped <- as.numeric(PHkap$snapped)

PHaep4$overall_mean_by_dist[PHaep4$overall_mean_by_dist == 0] <- NA
PHgno$overall_mean_by_dist[PHgno$overall_mean_by_dist == 0] <- NA
PHkap$overall_mean_by_dist[PHkap$overall_mean_by_dist == 0] <- NA


plot(est_lmr_overall, PHaep4$overall_mean_by_dist, log="xy")
plot(est_lmr_overall, PHgno$overall_mean_by_dist, log="xy")
plot(est_lmr_overall, PHkap$overall_mean_by_dist, log="xy")


del1 <- log10(PHaep4$overall_mean_by_dist)-log10(est_lmr_overall)
del2 <- log10(PHgno$overall_mean_by_dist)-log10(est_lmr_overall)
del3 <- log10(PHkap$overall_mean_by_dist)-log10(est_lmr_overall)
boxplot(del1)
boxplot(del2)
boxplot(del3)

SA <- PHaep4[! is.na(PHaep4$overall_mean_by_dist) & abs(del1) > 0.1,]
SG <- PHgno[ ! is.na(PHgno$overall_mean_by_dist ) & abs(del2) > 0.1,]
SK <- PHkap[ ! is.na(PHkap$overall_mean_by_dist ) & abs(del3) > 0.1,]


for(s in unique(SA$comid)) {
  tmp <- SA[SA$comid == s,]
  for(d in tmp$decade) {
     message("AEP4: resetting quantiles and overall mean to NA for ",s," and ",d,
             " because means don't align within 0.1 log10")
     for(i in 6:32) {
        PHaep4[PHaep4$comid == s & PHaep4$decade == d, i] <- NA
     }
     PHaep4$overall_mean_by_dist[PHaep4$comid == s & PHaep4$decade == d] <- NA
  }
}
for(s in unique(SG$comid)) {
  tmp <- SG[SG$comid == s,]
  for(d in tmp$decade) {
     message("GNO: resetting quantiles and overall mean to NA for ",s," and ",d,
             " because means don't align within 0.1 log10")
     for(i in 6:32) {
        PHgno[PHgno$comid == s & PHgno$decade == d, i] <- NA
     }
     PHgno$overall_mean_by_dist[PHgno$comid == s & PHgno$decade == d] <- NA
  }
}
for(s in unique(SK$comid)) {
  tmp <- SK[SK$comid == s,]
  for(d in tmp$decade) {
     message("KAP: resetting quantiles and overall mean to NA for ",s," and ",d,
             " because means don't align within 0.1 log10")
     for(i in 6:32) {
        PHkap[PHkap$comid == s & PHkap$decade == d, i] <- NA
     }
     PHkap$overall_mean_by_dist[PHkap$comid == s & PHkap$decade == d] <- NA
  }
}

plot(est_lmr_overall, PHaep4$overall_mean_by_dist, log="xy")
plot(est_lmr_overall, PHgno$overall_mean_by_dist,  log="xy")
plot(est_lmr_overall, PHkap$overall_mean_by_dist,  log="xy")

PHaep4[is.na(PHaep4$comid),]
PHgno[ is.na(PHgno$comid),]
PHkap[ is.na(PHkap$comid),]

as.data.frame(PHaep4[! is.na(PHaep4$f99.98) & PHaep4$f99.98 == 0,])
as.data.frame(PHgno[! is.na(PHgno$f99.98) & PHgno$f99.98 == 0,])
as.data.frame(PHkap[! is.na(PHkap$f99.98) & PHkap$f99.98 == 0,])


for(i in 1:length(PHaep4$comid)) {
  if(is.na(PHaep4$overall_mean_by_dist[i])) {
     for(j in 6:32) PHaep4[i,j] <- NA
  }
  if(is.na(PHgno$overall_mean_by_dist[i])) {
     for(j in 6:32) PHgno[i,j] <- NA
  }
  if(is.na(PHkap$overall_mean_by_dist[i])) {
     for(j in 6:32) PHkap[i,j] <- NA
  }
}

as.data.frame(PHaep4[! is.na(PHaep4$f99.98) & PHaep4$f99.98 == 0,])
as.data.frame(PHgno[! is.na(PHgno$f99.98) & PHgno$f99.98 == 0,])
as.data.frame(PHkap[! is.na(PHkap$f99.98) & PHkap$f99.98 == 0,])


write_feather(PHaep4, "all_huc_fdclmr_aep.feather")
write_feather(PHgno,  "all_huc_fdclmr_gno.feather")
write_feather(PHkap,  "all_huc_fdclmr_kap.feather")


PHaep4 <- read_feather("all_huc_fdclmr_aep.feather")
PHgno  <- read_feather("all_huc_fdclmr_gno.feather")
PHkap  <- read_feather("all_huc_fdclmr_kap.feather")


spPHaep4 <- SpatialPointsDataFrame(cbind(PHaep4$dec_long_va,
                                       PHaep4$dec_lat_va),
                        data=PHaep4, proj4string=LATLONG)
spPHaep4 <- spTransform(spPHaep4, ALBEA)
spPHgno <- SpatialPointsDataFrame(cbind(PHgno$dec_long_va,
                                       PHgno$dec_lat_va),
                        data=PHgno, proj4string=LATLONG)
spPHgno <- spTransform(spPHgno, ALBEA)
spPHkap <- SpatialPointsDataFrame(cbind(PHkap$dec_long_va,
                                       PHkap$dec_lat_va),
                        data=PHkap, proj4string=LATLONG)
spPHkap <- spTransform(spPHkap, ALBEA)

pdf("junk.pdf", useDingbats=FALSE)
  plot(spPHgno, lwd=0.2, pch=1, cex=0.5); mtext("AEP4 distribution")
  plot(spPHaep4[is.na(spPHaep4$overall_mean_by_dist),], add=TRUE, col=2, pch=2, cex=0.4)
  plot(spPHgno, lwd=0.2, pch=1, cex=0.5); mtext("GNO distribution")
  plot(spPHgno[is.na(spPHgno$overall_mean_by_dist),], add=TRUE, col=3, pch=3, cex=0.4)
  plot(spPHgno, lwd=0.2, pch=1, cex=0.5); mtext("KAP distribution")
  plot(spPHkap[is.na(spPHkap$overall_mean_by_dist),], add=TRUE, col=4, pch=4, cex=0.4)
dev.off()


