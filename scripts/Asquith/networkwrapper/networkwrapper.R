suppressWarnings(library(feather))
library(sp)
library(mgcv)
library(kernlab)

LATLONG <- paste0("+proj=longlat +ellps=GRS80 ",
                  "+datum=NAD83 +no_defs +towgs84=0,0,0")
LATLONG <- CRS(LATLONG)
ALBEA <- paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 ",
                "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
ALBEA <- CRS(ALBEA)

H   <- read.table("../../../scibase/huc12/all_huc12_covariates.csv", sep=",", header=TRUE,
                  colClasses="character")
length(H$comid)                           # [1] 55320 (20190311)

#SOL <- read.table("../../../scibase/huc12/all_huc12_solar.csv", sep=",", header=TRUE,
#                    colClasses="character")
#length(SOL$comid)                         # [1] 55128 (20190311)
# # ScienceBase is missing rows (55320-55128)/6 --- > 32 gages
SOL <- read_feather("../../../data/huc12/all_huc12_solar.feather")
H <- merge(H,SOL)

H$dec_lat_va  <-       as.numeric(H$dec_lat_va)
H$dec_long_va <-       as.numeric(H$dec_long_va)
H$basin_area  <- log10(as.numeric(H$basin_area))
H$basin_slope <- log10(as.numeric(H$basin_slope))
H$elev_mean   <- log10(as.numeric(H$elev_mean)+100)
H$ppt_mean    <- log10(as.numeric(H$ppt_mean))
H$decade      <- as.factor(H$decade)

GAG <- read.table("../../../scibase/gage/all_gage_covariates.csv", sep=",", header=TRUE,
                  colClasses="character")
STA <- read.table("../../../scibase/gage/all_gage_flow_stats.csv", sep=",", header=TRUE,
                  colClasses="character")
SOL <- read_feather("../../../data/gage/all_gage_solar.feather")

G <- merge(GAG,STA); rm(GAG); rm(STA)
G$f90 <- as.numeric(G$f90)
G$dec_lat_va  <-       as.numeric(G$dec_lat_va)
G$dec_long_va <-       as.numeric(G$dec_long_va)
G$basin_area  <- log10(as.numeric(G$basin_area))
G$basin_slope <- log10(as.numeric(G$basin_slope))
G$elev_mean   <- log10(as.numeric(G$elev_mean)+100)
G$ppt_mean    <- log10(as.numeric(G$ppt_mean))
G$decade <- as.factor(G$decade)
G$dni_ann <- NA
for(site in G$site_no) {
  G$dni_ann[G$site_no == site] <- rep(SOL$dni_ann[SOL$site_no == site][1],
                  length(G$dni_ann[G$site_no == site]))
}

H <- SpatialPointsDataFrame(cbind(H$dec_long_va, H$dec_lat_va),
                       data=H, proj4string=LATLONG)
H <- spTransform(H, ALBEA)

G <- SpatialPointsDataFrame(cbind(G$dec_long_va, G$dec_lat_va),
                       data=G, proj4string=LATLONG)
G <- spTransform(G, ALBEA)

XY <- coordinates(H)/1000
H$xkm <- XY[,1]; H$ykm <- XY[,2]; rm(XY)
XY <- coordinates(G)/1000
G$xkm <- XY[,1]; G$ykm <- XY[,2]; rm(XY)

Z <- log10(G$f90); tmp <- G[is.finite(Z),]; Z <- Z[is.finite(Z)]
tmp$Z <- Z; rm(Z)
#tmp <- tmp[tmp$xkm > -186839/1000,] # from locator()
GAM <-  gam(Z~s(basin_area)+s(basin_slope)+s(elev_mean)+s(ppt_mean)+s(xkm, ykm)+dni_ann+decade, data=tmp)
SVM <- ksvm(Z~  basin_area +  basin_slope +  elev_mean +  ppt_mean +  xkm+ ykm +dni_ann+decade, data=tmp)

h <- H[H$decade == 2000,]
h <- h[h$xkm > -186839/1000,] # from locator()
pGAM <- as.vector(predict(GAM, newdata=h))
pSVM <- as.vector(predict(SVM, newdata=h))

pdf("networktest_all.pdf", useDingbats=FALSE)
  h$DIFF <- abs(pGAM-pSVM)#/((pGAM+pSVM)/2)
  quas <- quantile(h$DIFF, probs=c(0.90,0.95,0.98))
  plot(H[H$decade == 2000,], col=NA)
  #plot(H[H$decade == 2000 & H$xkm <= -186839/1000,], add=TRUE, pch=1, col=grey(0.8), cex=H$basin_area/4, lwd=0.3)
  plot(H[H$decade == 2000,], add=TRUE, pch=1, col=grey(0.8), cex=H$basin_area/4, lwd=0.3)
  set <- (h$DIFF < quas[1])
  plot(h[set,], add=TRUE, pch=1, col=grey(0.4),  cex=h$basin_area/4, lwd=0.3)
  set <- (h$DIFF >= quas[1]) & (h$DIFF < quas[2])
  plot(h[set,], add=TRUE, pch=1, lwd=0.6, col=3, cex=h$basin_area/4)
  set <- (h$DIFF >= quas[2]) & (h$DIFF < quas[3])
  plot(h[set,], add=TRUE, pch=1, lwd=0.6, col=4, cex=h$basin_area/4)
  set <- (h$DIFF >= quas[3])
  plot(h[set,], add=TRUE, pch=1, lwd=0.6, col=2, cex=h$basin_area/4)
dev.off()







