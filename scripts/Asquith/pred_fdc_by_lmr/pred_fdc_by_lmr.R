library(sp)
library(feather)

LATLONG <- paste0("+proj=longlat +ellps=GRS80 ",
                  "+datum=NAD83 +no_defs +towgs84=0,0,0")
LATLONG <- sp::CRS(LATLONG)
ALBEA <- paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 ",
                "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
ALBEA <- sp::CRS(ALBEA)

COV <- read_feather("../../../data/huc12/all_huc12_covariates.feather")
SO <- read_feather("../../../data/huc12/all_huc12_solar.feather")
COV <- cbind(COV, SO)
spCOV <- SpatialPointsDataFrame(cbind(COV$dec_long_va,COV$dec_lat_va), data=COV,
                                proj4string=LATLONG)
spCOV <- spTransform(spCOV, ALBEA)
XY <- coordinates(spCOV)
spCOV$x <- spCOV$east <- XY[,1]/1000; spCOV$y <- spCOV$north <- XY[,2]/1000
rm(COV, XY)
spCOV$comid.1 <- NULL
spCOV$huc12.1 <- NULL
spCOV$dec_long_va.1 <- NULL
spCOV$dec_lat_va.1 <- NULL
spCOV$decade.1 <- NULL


library(lmomco)
library(feather)
gPPLO <- read_feather("../../../results/gage/gam/all_gage_looest_pplo.feather")
gL1   <- read_feather("../../../results/gage/gam/all_gage_looest_L1.feather")
gT2   <- read_feather("../../../results/gage/gam/all_gage_looest_T2.feather")
gT3   <- read_feather("../../../results/gage/gam/all_gage_looest_T3.feather")
gT4   <- read_feather("../../../results/gage/gam/all_gage_looest_T4.feather")
gT5   <- read_feather("../../../results/gage/gam/all_gage_looest_T5.feather")

hPPLO <- read_feather("../../../results/huc12/gam/all_huc12_pplo.feather")
hL1   <- read_feather("../../../results/huc12/gam/all_huc12_L1.feather")
hT2   <- read_feather("../../../results/huc12/gam/all_huc12_T2.feather")
hT3   <- read_feather("../../../results/huc12/gam/all_huc12_T3.feather")
hT4   <- read_feather("../../../results/huc12/gam/all_huc12_T4.feather")
hT5   <- read_feather("../../../results/huc12/gam/all_huc12_T5.feather")

G   <- read_feather("../../../data/gage/all_gage_data.feather")

FF <- 1/(3653+1); FF <- seq(qnorm(FF), qnorm(1-FF), 0.001)
FF <- pnorm(FF); FF <- c(0.00001, 0.99999, FF)
FF <- c(0.0001, 0.0015, 0.0002, FF) # for deep left tail smoothness
FF <- sort(unique(FF))
#huc12 <- sample(hL1$huc12, 1)
#comid <- sample(hL1$comid, 1)
site <- "08167000"
huc12 <- "121002010302" # at Guad Comfort
comid <- "3588922" # at Guad Comfort
plot(spCOV, lwd=0.2, col=8)
plot(spCOV[spCOV$huc12 == huc12,], col=2, add=TRUE)
plot(spCOV[spCOV$comid == comid,], col=4, add=TRUE)

COV <- spCOV@data
hmp <- COV[COV$huc12 == huc12,]
hmp_area <- hmp$basin_area[1]
cmp <- COV[COV$comid == comid,]
cmp_area <- cmp$basin_area[1]

gage <- G[G$site_no == site,]
gage_area <- gage$basin_area[1]
fdc <- as.data.frame(gage[,c(11:37)])
ff <- as.numeric(gsub("f","",names(fdc)))/100
ff <- qnorm(ff)
fdc <- as.data.frame(t(as.matrix(fdc)))
names(fdc) <- as.character(gage$decade)
phi <- 0.9 # Asquith and others (2006): Drainage-area ratio
duan <- gL1$bias_corr[1]

zero <- 0.001*(0.3048^3)
pdf("pred_fdc.pdf", useDingbats=FALSE, height=6, width=7)
opts <- par(); par(las=1, mgp=c(3,0.5,0))

plot(1,1, xlim=qnorm(c(0.0001, 0.9999)), ylim=c(0.001,500), xaxt="n", yaxt="n", xlab="",
          log="y", type="n", xaxs="i", yaxs="i",
          ylab="Projected and estimated mean streamflow, cms")
add.lmomco.axis(las=2, tcl=0.5, side.type="NPP", cex=0.8, case="lower")
add.log.axis(side=2, tcl=0.8*abs(par()$tcl), two.sided=TRUE)
#add.log.axis(logs=c(1),   tcl=-0.5*abs(par()$tcl), side=2, two.sided=TRUE)
add.log.axis(logs=c(1),   tcl=+1.3*abs(par()$tcl), side=2, two.sided=TRUE)
add.log.axis(logs=c(1, 2, 4, 6), side=2, make.labs=TRUE, las=1, label="")

ks <- c(0.6,0.8,1.0,1.2,1.4,1.6)-0.1
col <- .8; pch <- 3
for(i in 1:length(names(fdc))) {
  gmp <- fdc[,i]; #gmp[gmp == 0] <- 0.001
  lines(ff,  gmp*(cmp_area/gage_area)^phi, col=grey(col))
  points(ff, gmp*(cmp_area/gage_area)^phi, pch=21, bg=grey(col), lwd=0.7, cex=ks[i])
  col <- col - .1; #pch <- ifelse(pch == 3, 4, 3)
}

text(-3.4, 300, paste0("USGS streamgage ",site," ('next to' COMID ",comid," has area of ",gage_area," square kilometers"), pos=4, cex=0.7)
text(-3.4, 200, paste0("COMID ",comid," has area of ",cmp_area," square kilometers"), pos=4, cex=0.7)
text(-3.4, 60, paste0("'Projected' flow-duration curve (FDC) for COMID = \n   gaged FDC * ",
                       "(COMID basin_area / streamgage basin_area)^0.9 = flow,\n   in cubic meters per second (cms)"),
                pos=4, cex=0.7)
text(-1.0,.003, paste0("A Duan smearing estimator (Helsel and Hirsch, 2002) of ",
                     round(duan, digits=4), " was\n",
                     "used to correct retransformation bias in the mean (first L-moment)\n",
                     "prior to fitting of the distributions to the L-moments."),
             pos=4, cex=0.7)
txt <- c(paste0("Projected decadal FDC where circle size and greyness increases\n",
                "by decade from 1950s to 2000s"))
legend(-3.3, 40, txt, pt.bg=grey(0.7), bty="n", lwd=0.7, pch=21, cex=0.7, pt.cex=0.9)



for(decade in cmp$decade) {
  lwd <- ifelse(decade == "1950", 0.75,
         ifelse(decade == "1960", 1.15,
         ifelse(decade == "1970", 1.55,
         ifelse(decade == "1980", 0.75,
         ifelse(decade == "1990", 1.15,
         ifelse(decade == "2000", 1.55, 2))))))
  lty <- ifelse(decade == "1950", 2,
         ifelse(decade == "1960", 2,
         ifelse(decade == "1970", 2,
         ifelse(decade == "1980", 1,
         ifelse(decade == "1990", 1,
         ifelse(decade == "2000", 1, 3))))))
  message(decade, "  mu=", appendLF=FALSE)
  lmr <- vec2lmom(c(hL1$est_L1[hL1$comid == comid & hL1$decade == decade]*duan,
                    hT2$est_T2[hT2$comid == comid & hT2$decade == decade],
                    hT3$est_T3[hT3$comid == comid & hT3$decade == decade],
                    hT4$est_T4[hT4$comid == comid & hT4$decade == decade],
                    hT5$est_T5[hT5$comid == comid & hT5$decade == decade]
                    ), lscale=FALSE)

  aep4 <- paraep4(lmr, snap.tau4=FALSE)
  #print(aep4)
  aep4 <- paraep4(lmr, snap.tau4=TRUE)
  my.pplo <- hPPLO$est_pplo[hPPLO$comid == comid & hPPLO$decade == decade]
  nep1 <- qnorm(f2f(FF, pp=my.pplo))
  qua1 <- qua1o <- qlmomco(f2flo(FF, pp=my.pplo), aep4)
  qua1[qua1 <= zero] <- zero
  lines(c(min(nep1),nep1),c(zero,qua1), col="red", lwd=lwd, lty=lty)
  gno <- pargno(lmr)
  nep2 <- qnorm(f2f(FF, pp=my.pplo))
  qua2 <- qua2o <- qlmomco(f2flo(FF, pp=my.pplo), gno)
  qua2[qua2 <= zero] <- zero
  lines(c(min(nep2),nep2),c(zero,qua2), col="#006F41", lwd=lwd, lty=lty)
  kap <- parkap(lmr, snap.tau4=FALSE)
  if(kap$ifail > 0) message(decade, " decade above GLO for KAP")
  kap <- parkap(lmr, snap.tau4=TRUE)
  nep3 <- qnorm(f2f(FF, pp=my.pplo))
  qua3 <- qua3o <- qlmomco(f2flo(FF, pp=my.pplo), kap)
  qua3[qua3 <= zero] <- zero
  lines(c(min(nep3),nep3),c(zero,qua3), col="blue", lwd=lwd, lty=lty)
  df1 <- data.frame(nep=pnorm(nep1), qua1=qua1o)
  df2 <- data.frame(nep=pnorm(nep2), qua2=qua2o)
  df3 <- data.frame(nep=pnorm(nep3), qua3=qua3o)
  df <- merge(df1, merge(df2, df3, all=TRUE), all=TRUE)
  df$FF <- pnorm(unique(sort(c(nep1, nep2, nep3))))
  df$qua1[df$qua1 <= zero] <- zero
  df$qua2[df$qua2 <= zero] <- zero
  df$qua3[df$qua3 <= zero] <- zero
  df$Qfo <- df$Qfp <- 10^((log10(df$qua1) + log10(df$qua2) + log10(df$qua3))/3)
  #df$Qfo <- df$Qfp <- sapply(1:length(df$qua1), function(i) median(c(df$qua1[i], df$qua2[i], df$qua3[i])))
  df$Qfp[df$Qfp <= 0] <- 0
  row1 <- data.frame(nep=1/(3653+1), qua1=df$qua1[1],
                     qua2=df$qua2[1], qua3=df$qua3[1], FF=1/(3653+1),
                      Qfp=df$Qfp[1], Qfo=df$Qfo[1])
  df <- merge(df, row1, all=TRUE)
  #lines(qnorm(df$nep), df$Qfp, col=6, lwd=3, lty=1)

  comboMu <- function(f, combo) approx(combo$nep,combo$Qfo,xout=f,rule=2)$y
  mu <- integrate(comboMu, lower=0, upper=1, combo=df)$value
  message(mu)

  print(lmr$ratios[4:5])
  print(theoLmoms(aep4)$ratios[4:5])
  print(theoLmoms(gno)$ratios[4:5])
  print(theoLmoms(kap)$ratios[4:5])

}
mtext(paste0("COMID: ",comid))
txt <- c("Asymmetric exponential power distribution (four parameter)",
         "Generalized normal distribution (three parameter log-Normal)",
         "Kappa distribution (four parameter)")
legend(-3.4, 16, txt, col=c("red","#006F41","blue"), bty="n", lwd=1, cex=0.7)
text(-3.0, 1.9, paste0("Solid lines are 1980s to 2000s, dashed lines are 1950s to 1970s,\n",
                     "   and line width increases for the distributions."),
              cex=0.6, pos=4)
suppressWarnings(par(opts))
dev.off()


