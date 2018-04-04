library(feather)
library(mgcv)
library(sp)
library(GISTools)
library(RColorBrewer)

load("DEMO.RData")
load("Models.RData")
load(file.choose()) # GulfStates.RData
load(file.choose()) # spDNI_1998to2009.RData
load(file.choose()) # "spRESTORE_MGCV_BND.RData"
bnd <- list(x=bnd_poly_aea[,1]/1000, y=bnd_poly_aea[,2]/1000)
bnd <- list(bnd)

LATLONG <- paste0("+init=epsg:4269 +proj=longlat +ellps=GRS80 ",
                  "+datum=NAD83 +no_defs +towgs84=0,0,0")
LATLONG <- sp::CRS(LATLONG)
ALBEA <- paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 ",
                "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
ALBEA <- sp::CRS(ALBEA)

east_grids  <- seq(80,100,by=2)
north_grids <- seq(26,38, by=2)

gx <- gy <- vector(mode="numeric")
for(i in 1:length(north_grids)) {
   gy <- c(gy, rep(north_grids[i], length(east_grids)))
   gx <- c(gx, east_grids)
}
GL <- SpatialPoints(cbind(-gx, gy), proj4string=LATLONG)
GL <- spTransform(GL, ALBEA)
XY <- coordinates(GL)
x <- XY[,1]/1000; y <- XY[,2]/1000
#ind <- mgcv::inSide(bnd,x,y)
#XY <- XY[ind,]
GL <- SpatialPoints(cbind(x,y), proj4string=ALBEA)
ix <- 1:length(x)
plot(GL, pch=1, col=2)
text(XY[,1],XY[,2], ix)
GL <- GL[-c(1,3:9, 12, 14:20, 23, 34, 45)]

COV <- COVo <- read_feather(file.choose()) # "all_huc12_covariates.feather"
length(COVo$comid)                           # [1] 59070
length(COV$comid[is.na(COV$comid)])          # [1] 1458
COV <- COV[! is.na(COV$comid),]
length(COV$comid[is.na(COV$acc_basin_area)]) # [1] 330
COV <- COV[! is.na(COV$acc_basin_area), ]
length(COV$comid[is.na(COV$area_sqkm)])      # [1] 210

spCOV <- SpatialPointsDataFrame(cbind(COV$lon,COV$lat), data=COV,
                                proj4string=LATLONG)
spCOV <- spTransform(spCOV, ALBEA)
XY <- coordinates(spCOV)
spCOV$x <- spCOV$east <- XY[,1]/1000; spCOV$y <- spCOV$north <- XY[,2]/1000; rm(XY)



SO <- over(spCOV, spDNI_1998to2009)
spCOV$ANN_DNI <- SO$ANN_DNI
spCOV$JAN <- SO$JAN; spCOV$FEB <- SO$FEB
spCOV$MAR <- SO$MAR; spCOV$APR <- SO$APR
spCOV$MAY <- SO$MAY; spCOV$JUN <- SO$JUN
spCOV$JUL <- SO$JUL; spCOV$AUG <- SO$AUG
spCOV$SEP <- SO$SEP; spCOV$OCT <- SO$OCT
spCOV$NOV <- SO$NOV; spCOV$DEC <- SO$DEC
rm(SO)



COV$acc_nid_storage[COV$acc_basin_area == 0]
COV$acc_norm_storage[COV$acc_basin_area == 0]


# storages are in acre-ft
# 1 km2 = 247.104393047 acres
spCOV$flood_storage <- spCOV$acc_nid_storage - spCOV$acc_norm_storage
spCOV$flood_storage <- spCOV$flood_storage/(spCOV$acc_basin_area*247.104393047)
spCOV$flood_storage[is.na(spCOV$flood_storage)] <- 0
spCOV[spCOV$flood_storage < 0,]; # plot(spCOV); plot(spCOV[spCOV$flood_storage < 0,], add=TRUE, col=2)
spCOV$flood_storage <- abs(spCOV$flood_storage)
spCOV$flood_storage <- log10(spCOV$flood_storage+.01)
plot(qnorm(pp(spCOV$flood_storage)), sort(spCOV$flood_storage), type="l")
#lines(qnorm(pp(DD$flood_storage)), sort(DD$flood_storage), col=2)
length(spCOV$comid[spCOV$flood_storage >= 0.52]) # This is about 3.3 feet of watershed depth storage
# [1] 35
#summary(spCOV$flood_storage[spCOV$flood_storage >= 0.52])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.5228  0.5546  0.7193  0.7827  0.8425  1.5351

spCOV <- spCOV[spCOV$flood_storage < 0.52,]


spCOV$acc_basin_slope[spCOV$acc_basin_slope == 0] <- 0.0005

length(spCOV$comid[spCOV$acc_basin_area == 0])
spCOV <- spCOV[spCOV$acc_basin_area != 0,]

spCOV$ppt_mean        <- log10(spCOV$ppt_mean)
spCOV$temp_mean       <- log10(spCOV$temp_mean)
spCOV$acc_basin_area  <- log10(spCOV$acc_basin_area)
spCOV$acc_basin_slope <- log10(spCOV$acc_basin_slope)
length(spCOV$comid[! is.finite(spCOV$ppt_mean)])
length(spCOV$comid[! is.finite(spCOV$temp_mean)])
length(spCOV$comid[! is.finite(spCOV$acc_basin_area)])
length(spCOV$comid[! is.finite(spCOV$acc_basin_slope)])



spCOV$decade <- as.factor(spCOV$decade)
spCOV$bedperm <- as.factor(spCOV$bedperm)
spCOV$aquifers <- as.factor(spCOV$aquifers)
spCOV$soller <- as.factor(spCOV$soller)
spCOV$hlr <- as.factor(spCOV$hlr)
spCOV$ecol3 <- as.factor(spCOV$ecol3)
spCOV$physio <- as.factor(spCOV$physio)
spCOV$statsgo <- as.factor(spCOV$statsgo)

spCOV$barren <- 2*asin(sqrt(spCOV$barren/100))
spCOV$cultivated_cropland <- 2*asin(sqrt(spCOV$cultivated_cropland/100))
spCOV$deciduous_forest    <- 2*asin(sqrt(spCOV$deciduous_forest/100))
spCOV$developed           <- 2*asin(sqrt(spCOV$developed/100))
spCOV$evergreen_forest    <- 2*asin(sqrt(spCOV$evergreen_forest/100))
spCOV$grassland           <- 2*asin(sqrt(spCOV$grassland/100))
spCOV$hay_pasture         <- 2*asin(sqrt(spCOV$hay_pasture/100))
spCOV$herbaceous_wetland  <- 2*asin(sqrt(spCOV$herbaceous_wetland/100))
spCOV$mixed_forest        <- 2*asin(sqrt(spCOV$mixed_forest/100))
# Note perennial_ice_snow is 0.00 throughout
spCOV$perennial_ice_snow  <- 2*asin(sqrt(spCOV$perennial_ice_snow/100))
spCOV$shrubland           <- 2*asin(sqrt(spCOV$shrubland/100))
spCOV$water               <- 2*asin(sqrt(spCOV$water/100))
spCOV$woody_wetland       <- 2*asin(sqrt(spCOV$woody_wetland/100))

message("REMOVING nodata (Bed Permeability)")
length(spCOV$comid[spCOV$bedperm == "nodata"]) # [1] 198
spCOV <- spCOV[spCOV$bedperm != "nodata",]

message("REMOVING ecol3_37, ecol3_72, nodata (Ecoregion)")
length(spCOV$comid[spCOV$ecol3 == "ecol3_37"]) # [1] 198
length(spCOV$comid[spCOV$ecol3 == "ecol3_72"]) # [1] 12
length(spCOV$comid[spCOV$ecol3 == "nodata"])   # [1] 30



map_annotation <- function() {
  ss <- list(x=-90000, y=320000)
  SpatialPolygonsRescale(layout.scale.bar(height=.07), offset=ss, scale=200*1000,
                         fill=c("transparent", "black"), plot.grid=FALSE, lwd=0.4)
  ss <- list(x=-90000, y=303000)
  SpatialPolygonsRescale(layout.scale.bar(height=.08), offset=ss, scale=100*1609.344,
                         fill=c("transparent", "black"), plot.grid=FALSE, lwd=0.4)
  xx <- -122000
  sl <- list(x=xx, y=350000)
  sr <- list(x=xx+200*1000, y=350000)
  text(sl, "0",xx, cex=0.6, pos=4); text(sr, "200 kilometers", cex=0.5, pos=4)
  sl <- list(x=xx, y=283000)
  sr <- list(x=xx+100*1609.344, y=283000)
  text(sl, "0", cex=0.6, pos=4); text(sr, "100 miles", cex=0.5, pos=4)

  txt <- paste0("Albers Equal Area Projection\n",
                "North American Datum of 1983\n",
                "Base modified from U.S. Geological Survey digital data, 1:24,000")
  text(-420000, 1590000, txt, pos=4, cex=0.45)
  plot(GulfStates_modified, add=TRUE, lwd=.4, lty=2)
  STATES <- c("Texas", "Oklahoma", "Missouri", "Arkansas", "Louisiana", "Mississippi",
              "Tennessee", "Kentucky", "Alabama", "Georgia", "Florida")
  STATES <- data.frame(easting=c(-440000, -202900.4,  178000,  178000,  279961.5,
                                 430000,  710000,  710000,  710000,
                                 1100000, 1270000),
                       northing=c(955139.0, 1400000, 1558716.4, 1400000,  795368.5,
                                  1170000, 1450000, 1600000, 1325000,
                                  1170000, 800000),
                       state=STATES)
  text(STATES$easting, STATES$northing, STATES$state, pos=4, cex=0.7, col=grey(0.3))
  plot(GL, lwd=0.4, col=grey(0.22), add=TRUE)
}

map_base <- function() {
  par(lend=1, ljoin=1)
  plot(spCOV, pch=NA); plot(GulfStates, add=TRUE, lty=0, col=grey(0.95))
  polygon(bnd[[1]]$x*1000,bnd[[1]]$y*1000, col=grey(1), lwd=.7)
}

map_sebase <- function() {
  k <- 0; ks <- c(0.6,0.8,1.0,1.2,1.4,1.6)
  for(d in sort(unique(D$decade))) {
    k <- k + 1
    plot(D[D$decade == d,], pch=1, lwd=0.7, cex=ks[k], col=grey(0.72), add=TRUE)
  }
}

choropleth_decade <-
  function(data, cuts, x=NA, rev=FALSE, trans=function(t) {t}) {
    env <- as.environment(as.list(slot(data, "data")))
    if(x == "T2") {
      y <- get("L2", envir=env)
      x <- get("L1", envir=env)
      x <- y/x
    } else {
      x <- get(x, envir=env)
    }
    decade <- get("decade", envir=env)
    cols <- add.alpha(brewer.pal(10,"Spectral"),.7)
    if(rev) cols <- rev(cols)
    k <- 0; ks <- c(0.6,0.8,1.0,1.2,1.4,1.6)
    for(d in sort(unique(decade))) {
      k <- k + 1
      tmp <- trans(x[decade == d])
      shades <- auto.shading(tmp, cutter=cuts, n=9, cols=cols)
      choropleth(data[data$decade == d,], tmp, pch=1, lwd=0.7, cex=ks[k], shading=shades, add=TRUE)
    }
  }

choropleth_cov <-
  function(data, cuts, x=NA, decade="2000", rev=FALSE, trans=function(t) {t}) {
    data <- data[data$decade == decade,]
    env <- as.environment(as.list(slot(data, "data")))
    x <- get(x,   envir=env)
    cols <- add.alpha(brewer.pal(10,"Spectral"),.7)
    if(rev) cols <- rev(cols)
    tmp <- trans(x)
    shades <- auto.shading(tmp, cutter=cuts, n=9, cols=cols)
    choropleth(data, tmp, pch=16, cex=0.4, shading=shades, add=TRUE)
    return(shades)
  }

legend_est <- function(gage="", title="", note=TRUE, shades=NA, itgage=TRUE, ...) {
  tx1 <- paste0("U.S. Geological Survey streamflow-gaging station:\n",
                "Symbol colored* by observed decadal ",gage)
  tx2 <- paste0("'COMID' location in the National Hydrography Dataset\n",
                "version 2: Symbol colored by estimated decadal\n",gage)
  if(itgage) {
    legend(-80000, 595000, c(tx1, tx2), bty="n", cex=0.7, pt.cex=c(0.8,0.6),
           lwd=c(0.7,1), lty=c(0,0), pch=c(1,16), col=c("#3288BDE6","#D53E4FE6"), ...)
  } else {
    legend(-80000, 520000, c(tx2), bty="n", cex=0.7, pt.cex=c(0.6),
           lwd=c(1), lty=c(0), pch=c(16), col=c("#D53E4FE6"), ...)
  }
  choro.legend(850000, 690000, shades, cex=0.7, bty="n", box.col=grey(1), bg=grey(1),
               fmt="%g", xjust=0, title=title)
  if(note) {
    text(270000, 320000,
         "* Note that symbol size represents decade\n (1950 [smallest circle] through 2000 [largest circle]).",
         cex=.6, pos=4)
  }
}

setxt1 <- "standard error of fit"
setxt2 <- "Standard error of fit"










EPo <- predict(PPLO, newdata=spCOV, se.fit=TRUE)
EPo$se.fit <- EPo$se.fit*sqrt(PPLO$sig2)
EPo$fit[is.na(EPo$fit)] <- mean(EPo$fit, na.rm=TRUE)
EPo$fit[EPo$fit > log10(3653)] <- log10(3653)
EPo$se.fit[is.na(EPo$se.fit)] <- mean(EPo$se.fit, na.rm=TRUE)

spCOV$est_pplo <- (3653-10^EPo$fit)/3653
spCOV$est_pplo[is.na(spCOV$est_pplo)] <- 0
spCOV$est_pplo[spCOV$est_pplo < 1/3653] <- 0
spCOV$se.fit_est_pplo <- EPo$se.fit

quantile(spCOV$est_pplo, probs=(1:9)/10, na.rm=TRUE)
quantile(EPo$se.fit,    probs=(1:9)/10, na.rm=TRUE)

pploCuts <- function(x, n=9, ...) {
   labs <- 1:n
   cuts <- c(0, 0.005, .01, 0.04, .05, 0.1, 0.2, 0.4, 0.6)
   cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}

pploCutsSE <- function(x, n=9, ...) {
   labs <- 1:n
   cuts <- c(0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.019, 0.020, 0.024)
   cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}

EL1 <- predict(L1, newdata=spCOV, se.fit=TRUE)
#EL1$se.fit <- sqrt(L1$sig2) * sqrt((EL1$se.fit/sqrt(L1$sig2))^2+1)
EL1$fit[is.na(EL1$fit)] <- mean(EL1$fit, na.rm=TRUE)
EL1$se.fit[is.na(EL1$se.fit)] <- mean(EL1$se.fit, na.rm=TRUE)
spCOV$est_L1 <- EL1$fit
spCOV$se.fit_est_L1 <- EL1$se.fit

quantile(spCOV$est_L1, probs=(1:9)/10, na.rm=TRUE)
quantile(EL1$se.fit,    probs=(1:9)/10, na.rm=TRUE)

L1Cuts <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(1.2, 1.5, 1.6, 1.8, 1.9, 2.1, 2.3, 2.6, 3.1)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}

L1CutsSE <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.021)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}




ET2 <- predict(T2, newdata=spCOV, se.fit=TRUE)
ET2$fit[is.na(ET2$fit)] <- mean(ET2$fit, na.rm=TRUE)
ET2$se.fit[is.na(ET2$se.fit)] <- mean(ET2$se.fit, na.rm=TRUE)
spCOV$est_T2 <- ET2$fit
spCOV$se.fit_est_T2 <- ET2$se.fit

quantile(spCOV$est_T2, probs=(1:9)/10, na.rm=TRUE)
quantile(ET2$se.fit,    probs=(1:9)/10, na.rm=TRUE)

T2Cuts <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.55, 0.61, 0.65, 0.68, 0.71, 0.74, 0.77, 0.81, 0.86)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}

T2CutsSE <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.0080, 0.0085, 0.0089, 0.0093, 0.0098, 0.0103, 0.0111, 0.0124, 0.0140)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}


ET3 <- predict(T3, newdata=spCOV, se.fit=TRUE)
ET3$fit[is.na(ET3$fit)] <- mean(ET3$fit, na.rm=TRUE)
ET3$se.fit[is.na(ET3$se.fit)] <- mean(ET3$se.fit, na.rm=TRUE)
spCOV$est_T3 <- ET3$fit
spCOV$se.fit_est_T3 <- ET3$se.fit

quantile(spCOV$est_T3, probs=(1:9)/10, na.rm=TRUE)
quantile(ET3$se.fit,    probs=(1:9)/10, na.rm=TRUE)

T3Cuts <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.53, 0.58, 0.61, 0.64, 0.67, 0.70, 0.73, 0.77, 0.84)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}

T3CutsSE <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.0077, 0.0822, 0.0086, 0.0900, 0.0945, 0.0100, 0.0108, 0.0120, 0.0140)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}


ET4 <- predict(T4, newdata=spCOV, se.fit=TRUE)
ET4$fit[is.na(ET4$fit)] <- mean(ET4$fit, na.rm=TRUE)
ET4$se.fit[is.na(ET4$se.fit)] <- mean(ET4$se.fit, na.rm=TRUE)
spCOV$est_T4 <- ET4$fit
spCOV$se.fit_est_T4 <- ET4$se.fit

quantile(spCOV$est_T4, probs=(1:9)/10, na.rm=TRUE)
quantile(ET4$se.fit,    probs=(1:9)/10, na.rm=TRUE)

T4Cuts <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.34, 0.40, 0.44, 0.46, 0.50, 0.52, 0.56, 0.62, 0.72)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}

T4CutsSE <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.0081, 0.0087, 0.0091, 0.0095, 0.01000, 0.0105, 0.0114, 0.0126, 0.0145)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}


ET5 <- predict(T5, newdata=spCOV, se.fit=TRUE)
ET5$fit[is.na(ET5$fit)] <- mean(ET5$fit, na.rm=TRUE)
ET5$se.fit[is.na(ET5$se.fit)] <- mean(ET5$se.fit, na.rm=TRUE)
spCOV$est_T5 <- ET5$fit
spCOV$se.fit_est_T5 <- ET5$se.fit

quantile(spCOV$est_T5, probs=(1:9)/10, na.rm=TRUE)
quantile(ET5$se.fit,    probs=(1:9)/10, na.rm=TRUE)

T5Cuts <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.23, 0.229, 0.32, 0.36, 0.38, 0.41, 0.44, 0.50, 0.60)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}

T5CutsSE <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.0080, 0.0085, 0.0090, 0.0094, 0.0098,
            0.0104, 0.0112, 0.0124, 0.0143)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}


ET6 <- predict(T6, newdata=spCOV, se.fit=TRUE)
ET6$fit[is.na(ET6$fit)] <- mean(ET6$fit, na.rm=TRUE)
ET6$se.fit[is.na(ET6$se.fit)] <- mean(ET6$se.fit, na.rm=TRUE)
spCOV$est_T6 <- ET6$fit
spCOV$se.fit_est_T6 <- ET6$se.fit

quantile(spCOV$est_T6, probs=(1:9)/10, na.rm=TRUE)
quantile(ET6$se.fit,    probs=(1:9)/10, na.rm=TRUE)

T6Cuts <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.17, 0.22, 0.26, 0.29, 0.31, 0.33, 0.36, 0.41, 0.52)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}

T6CutsSE <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.0074, 0.0079, 0.0083, 0.0086, 0.0090, 0.0096, 0.0103, 0.0114, 0.0131)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}

EQ50 <- predict(Q50, newdata=spCOV, se.fit=TRUE)
EQ50$fit[is.na(EQ50$fit)] <- mean(EQ50$fit, na.rm=TRUE)
EQ50$se.fit[is.na(EQ50$se.fit)] <- mean(EQ50$se.fit, na.rm=TRUE)

spCOV$est_f50 <- EQ50$fit
spCOV$se.fit_est_f50 <- EQ50$se.fit

quantile(spCOV$est_f50,  probs=(1:9)/10, na.rm=TRUE)
quantile(EQ50$se.fit,    probs=(1:9)/10, na.rm=TRUE)

f50Cuts <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.21, 0.67, 0.94, 1.13, 1.31, 1.52, 1.76, 2.10, 2.65)
  cuts <- cuts[labs]
  names(cuts) <- paste("#", labs, sep = "")
  cuts
}

f50CutsSE <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.026, 0.028, 0.029, 0.031, 0.032, 0.034, 0.037, 0.041, 0.048)
  cuts <- cuts[labs]
  names(cuts) <- paste("#", labs, sep = "")
  cuts
}



EQ90 <- predict(Q90, newdata=spCOV, se.fit=TRUE)
EQ90$fit[is.na(EQ90$fit)] <- mean(EQ90$fit, na.rm=TRUE)
EQ90$se.fit[is.na(EQ90$se.fit)] <- mean(EQ90$se.fit, na.rm=TRUE)

spCOV$est_f90 <- EQ90$fit
spCOV$se.fit_est_f90 <- EQ90$se.fit

quantile(spCOV$est_f90,  probs=(1:9)/10, na.rm=TRUE)
quantile(EQ90$se.fit,    probs=(1:9)/10, na.rm=TRUE)

f90Cuts <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(1.3, 1.6, 1.8, 2.0, 2.2, 2.3, 2.6, 3.0, 3.5)
  cuts <- cuts[labs]
  names(cuts) <- paste("#", labs, sep = "")
  cuts
}

f90CutsSE <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.021, 0.022, 0.026)
  cuts <- cuts[labs]
  names(cuts) <- paste("#", labs, sep = "")
  cuts
}


EQ95 <- predict(Q95, newdata=spCOV, se.fit=TRUE)
EQ95$fit[is.na(EQ95$fit)] <- mean(EQ95$fit, na.rm=TRUE)
EQ95$se.fit[is.na(EQ95$se.fit)] <- mean(EQ95$se.fit, na.rm=TRUE)

spCOV$est_f95 <- EQ95$fit
spCOV$se.fit_est_f95 <- EQ95$se.fit

quantile(spCOV$est_f95,  probs=(1:9)/10, na.rm=TRUE)
quantile(EQ95$se.fit,    probs=(1:9)/10, na.rm=TRUE)

f95Cuts <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(1.6, 2.0, 2.1, 2.3, 2.4, 2.6, 2.8, 3.2, 3.7)
  cuts <- cuts[labs]
  names(cuts) <- paste("#", labs, sep = "")
  cuts
}

f95CutsSE <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.020, 0.022)
  cuts <- cuts[labs]
  names(cuts) <- paste("#", labs, sep = "")
  cuts
}

EQ98 <- predict(Q98, newdata=spCOV, se.fit=TRUE)
EQ98$fit[is.na(EQ98$fit)] <- mean(EQ98$fit, na.rm=TRUE)
EQ98$se.fit[is.na(EQ98$se.fit)] <- mean(EQ98$se.fit, na.rm=TRUE)

spCOV$est_f98 <- EQ98$fit
spCOV$se.fit_est_f98 <- EQ98$se.fit

quantile(spCOV$est_f98,  probs=(1:9)/10, na.rm=TRUE)
quantile(EQ98$se.fit,    probs=(1:9)/10, na.rm=TRUE)

f98Cuts <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(2.0, 2.3, 2.5, 2.6, 2.7, 2.9, 3.1, 3.4, 3.9)
  cuts <- cuts[labs]
  names(cuts) <- paste("#", labs, sep = "")
  cuts
}

f98CutsSE <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.020, 0.022)
  cuts <- cuts[labs]
  names(cuts) <- paste("#", labs, sep = "")
  cuts
}


EQ99 <- predict(Q99, newdata=spCOV, se.fit=TRUE)
EQ99$fit[is.na(EQ99$fit)] <- mean(EQ99$fit, na.rm=TRUE)
EQ99$se.fit[is.na(EQ99$se.fit)] <- mean(EQ99$se.fit, na.rm=TRUE)

spCOV$est_f99 <- EQ99$fit
spCOV$se.fit_est_f99 <- EQ99$se.fit

quantile(spCOV$est_f99,  probs=(1:9)/10, na.rm=TRUE)
quantile(EQ99$se.fit,    probs=(1:9)/10, na.rm=TRUE)

f99Cuts <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(2.3, 2.5, 2.7, 2.8, 2.9, 3.1, 3.3, 3.6, 4.1)
  cuts <- cuts[labs]
  names(cuts) <- paste("#", labs, sep = "")
  cuts
}

f99CutsSE <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.020, 0.021, 0.024)
  cuts <- cuts[labs]
  names(cuts) <- paste("#", labs, sep = "")
  cuts
}


EQ99p9 <- predict(Q99p9, newdata=spCOV, se.fit=TRUE)
EQ99p9$fit[is.na(EQ99p9$fit)] <- mean(EQ99p9$fit, na.rm=TRUE)
EQ99p9$se.fit[is.na(EQ99p9$se.fit)] <- mean(EQ99p9$se.fit, na.rm=TRUE)

spCOV$est_f99.9 <- EQ99p9$fit
spCOV$se.fit_est_f99.9 <- EQ99p9$se.fit

quantile(spCOV$est_f99.9,  probs=(1:9)/10, na.rm=TRUE)
quantile(EQ99p9$se.fit,    probs=(1:9)/10, na.rm=TRUE)

f99p9Cuts <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(2.9, 3.1, 3.2, 3.3, 3.4, 3.6, 3.8, 4.0, 4.4)
  cuts <- cuts[labs]
  names(cuts) <- paste("#", labs, sep = "")
  cuts
}

f99p9CutsSE <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.019, 0.021, 0.022, 0.023, 0.024, 0.025, 0.026, 0.029, 0.033)
  cuts <- cuts[labs]
  names(cuts) <- paste("#", labs, sep = "")
  cuts
}



#-----------------------------------------------------------------------
pdf("PPLOfit.pdf", useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    map_base()
    choropleth_decade(D, x="pplo", cuts=pploCuts, rev=TRUE)
    shades <- choropleth_cov(spCOV, decade=d, x="est_pplo", cuts=pploCuts, rev=TRUE)
    legend_est(gage="no flow fraction", title=paste0(d," decade\n","no flow fraction"),
               note=TRUE, shades=shades)
    map_annotation()
  }
  for(d in sort(unique(D$decade))) {
    map_base(); map_sebase()
    shades <- choropleth_cov(spCOV, decade=d, x="se.fit_est_pplo", cuts=pploCutsSE, rev=TRUE)
    legend_est(gage=setxt1, title=paste0(d," decade\n",setxt1), note=FALSE, shades=shades, itgage=FALSE)
    map_annotation()
  }
dev.off()
#-----------------------------------------------------------------------
pdf("L1fit.pdf", useDingbats=FALSE, width=11, height=10)
for(d in sort(unique(D$decade))) {
  map_base()
  choropleth_decade(D, x="L1", cuts=L1Cuts, trans=log10)
  shades <- choropleth_cov(spCOV, decade=d, x="est_L1", cuts=L1Cuts)
  legend_est(gage="mean streamflow", title=paste0(d," decade\n","mean streamflow, in log10(cfs)"),
             note=TRUE, shades=shades)
  map_annotation()
}
for(d in sort(unique(D$decade))) {
  map_base(); map_sebase()
  shades <- choropleth_cov(spCOV, decade=d, x="se.fit_est_L1", cuts=L1CutsSE, rev=TRUE)
  legend_est(gage=setxt1, title=paste0(d," decade\n",setxt1), note=FALSE, shades=shades, itgage=FALSE)
  map_annotation()
}
dev.off()
#-----------------------------------------------------------------------
pdf("T2fit.pdf", useDingbats=FALSE, width=11, height=10)
for(d in sort(unique(D$decade))) {
  map_base(); cuts <-
  choropleth_decade(D, x="T2", cuts=T2Cuts)
  shades <- choropleth_cov(spCOV, decade=d, x="est_T2", cuts=T2Cuts)
  legend_est(gage="L-CV of streamflow", title=paste0(d," decade\n","L-CV of streamflow"),
             note=TRUE, shades=shades)
  map_annotation()
}
for(d in sort(unique(D$decade))) {
  map_base(); map_sebase()
  shades <- choropleth_cov(spCOV, decade=d, x="se.fit_est_T2", cuts=T2CutsSE, rev=TRUE)
  legend_est(gage=setxt1, title=paste0(d," decade\n",setxt1), note=FALSE, shades=shades, itgage=FALSE)
  map_annotation()
}
dev.off()
#-----------------------------------------------------------------------
pdf("T3fit.pdf", useDingbats=FALSE, width=11, height=10)
for(d in sort(unique(D$decade))) {
  map_base()
  choropleth_decade(D, x="T3", cuts=T3Cuts)
  shades <- choropleth_cov(spCOV, decade=d, x="est_T3", cuts=T3Cuts)
  legend_est(gage="L-skew of streamflow", title=paste0(d," decade\n","L-skew of streamflow"),
             note=TRUE, shades=shades)
  map_annotation()
}
for(d in sort(unique(D$decade))) {
  map_base(); map_sebase()
  shades <- choropleth_cov(spCOV, decade=d, x="se.fit_est_T3", cuts=T3CutsSE, rev=TRUE)
  legend_est(gage=setxt1, title=paste0(d," decade\n",setxt1), note=FALSE, shades=shades, itgage=FALSE)
  map_annotation()
}
dev.off()
#-----------------------------------------------------------------------
pdf("T4fit.pdf", useDingbats=FALSE, width=11, height=10)
for(d in sort(unique(D$decade))) {
  map_base()
  choropleth_decade(D, x="T4", cuts=T4Cuts)
  shades <- choropleth_cov(spCOV, decade=d, x="est_T4", cuts=T4Cuts)
  legend_est(gage="L-kurtosis of streamflow", title=paste0(d," decade\n","L-kurtosis of streamflow"),
             note=TRUE, shades=shades)
  map_annotation()
}
for(d in sort(unique(D$decade))) {
  map_base(); map_sebase()
  shades <- choropleth_cov(spCOV, decade=d, x="se.fit_est_T4", cuts=T4CutsSE, rev=TRUE)
  legend_est(gage=setxt1, title=paste0(d," decade\n",setxt1), note=FALSE, shades=shades, itgage=FALSE)
  map_annotation()
}
dev.off()
#-----------------------------------------------------------------------
pdf("T5fit.pdf", useDingbats=FALSE, width=11, height=10)
for(d in sort(unique(D$decade))) {
  map_base()
  choropleth_decade(D, x="T5", cuts=T5Cuts)
  shades <- choropleth_cov(spCOV, decade=d, x="est_T5", cuts=T5Cuts)
  legend_est(gage="Tau5 of streamflow", title=paste0(d," decade\n","Tau5 of streamflow"),
             note=TRUE, shades=shades)
  map_annotation()
}
for(d in sort(unique(D$decade))) {
  map_base(); map_sebase()
  shades <- choropleth_cov(spCOV, decade=d, x="se.fit_est_T5", cuts=T5CutsSE, rev=TRUE)
  legend_est(gage=setxt1, title=paste0(d," decade\n",setxt1), note=FALSE, shades=shades, itgage=FALSE)
  map_annotation()
}
dev.off()
#-----------------------------------------------------------------------
pdf("T6fit.pdf", useDingbats=FALSE, width=11, height=10)
for(d in sort(unique(D$decade))) {
  map_base()
  choropleth_decade(D, x="T6", cuts=T6Cuts)
  shades <- choropleth_cov(spCOV, decade=d, x="est_T6", cuts=T6Cuts)
  legend_est(gage="T6 of streamflow", title=paste0(d," decade\n","Tau6 of streamflow"),
             note=TRUE, shades=shades)
  map_annotation()
}
for(d in sort(unique(D$decade))) {
  map_base(); map_sebase()
  shades <- choropleth_cov(spCOV, decade=d, x="se.fit_est_T6", cuts=T6CutsSE, rev=TRUE)
  legend_est(gage=setxt1, title=paste0(d," decade\n",setxt1), note=FALSE, shades=shades, itgage=FALSE)
  map_annotation()
}
dev.off()
#-----------------------------------------------------------------------
pdf("Q50fit.pdf", useDingbats=FALSE, width=11, height=10)
for(d in sort(unique(D$decade))) {
  map_base()
  choropleth_decade(D, x="f50", cuts=f50Cuts, trans=function(t) { log10(t+1) } )
  shades <- choropleth_cov(spCOV, decade=d, x="est_f50", cuts=f50Cuts)
  legend_est(gage="log10(Q50+1) of streamflow", title=paste0(d," decade\n","log10(Q50+1) of streamflow"),
             note=TRUE, shades=shades)
  map_annotation()
}
for(d in sort(unique(D$decade))) {
  map_base(); map_sebase()
  shades <- choropleth_cov(spCOV, decade=d, x="se.fit_est_f50", cuts=f50CutsSE, rev=TRUE)
  legend_est(gage=setxt1, title=paste0(d," decade\n",setxt1), note=FALSE, shades=shades, itgage=FALSE)
  map_annotation()
}
dev.off()
#-----------------------------------------------------------------------
pdf("Q90fit.pdf", useDingbats=FALSE, width=11, height=10)
for(d in sort(unique(D$decade))) {
  map_base()
  choropleth_decade(D, x="f90", cuts=f90Cuts, trans=function(t) { log10(t+1) } )
  shades <- choropleth_cov(spCOV, decade=d, x="est_f90", cuts=f90Cuts)
  legend_est(gage="log10(Q90+1) of streamflow", title=paste0(d," decade\n","log10(Q90+1) of streamflow"),
             note=TRUE, shades=shades)
  map_annotation()
}
for(d in sort(unique(D$decade))) {
  map_base(); map_sebase()
  shades <- choropleth_cov(spCOV, decade=d, x="se.fit_est_f90", cuts=f90CutsSE, rev=TRUE)
  legend_est(gage=setxt1, title=paste0(d," decade\n",setxt1), note=FALSE, shades=shades, itgage=FALSE)
  map_annotation()
}
dev.off()
#-----------------------------------------------------------------------
pdf("Q95fit.pdf", useDingbats=FALSE, width=11, height=10)
for(d in sort(unique(D$decade))) {
  map_base()
  choropleth_decade(D, x="f95", cuts=f95Cuts, trans=function(t) { log10(t+1) } )
  shades <- choropleth_cov(spCOV, decade=d, x="est_f95", cuts=f95Cuts)
  legend_est(gage="log(Q95+1) of streamflow", title=paste0(d," decade\n","log10(Q95+1) of streamflow"),
             note=TRUE, shades=shades)
  map_annotation()
}
for(d in sort(unique(D$decade))) {
  map_base(); map_sebase()
  shades <- choropleth_cov(spCOV, decade=d, x="se.fit_est_f95", cuts=f95CutsSE, rev=TRUE)
  legend_est(gage=setxt1, title=paste0(d," decade\n",setxt1), note=FALSE, shades=shades, itgage=FALSE)
  map_annotation()
}
dev.off()
#-----------------------------------------------------------------------
pdf("Q98fit.pdf", useDingbats=FALSE, width=11, height=10)
for(d in sort(unique(D$decade))) {
  map_base()
  choropleth_decade(D, x="f98", cuts=f98Cuts, trans=function(t) { log10(t+1) } )
  shades <- choropleth_cov(spCOV, decade=d, x="est_f98", cuts=f98Cuts)
  legend_est(gage="log(Q98+1) of streamflow", title=paste0(d," decade\n","log10(Q98+1) of streamflow"),
             note=TRUE, shades=shades)
  map_annotation()
}
for(d in sort(unique(D$decade))) {
  map_base(); map_sebase()
  shades <- choropleth_cov(spCOV, decade=d, x="se.fit_est_f98", cuts=f98CutsSE, rev=TRUE)
  legend_est(gage=setxt1, title=paste0(d," decade\n",setxt1), note=FALSE, shades=shades, itgage=FALSE)
  map_annotation()
}
dev.off()
#-----------------------------------------------------------------------
pdf("Q99fit.pdf", useDingbats=FALSE, width=11, height=10)
for(d in sort(unique(D$decade))) {
  map_base()
  choropleth_decade(D, x="f99", cuts=f99Cuts, trans=function(t) { log10(t+1) } )
  shades <- choropleth_cov(spCOV, decade=d, x="est_f99", cuts=f99Cuts)
  legend_est(gage="log(Q99+1) of streamflow", title=paste0(d," decade\n","log10(Q99+1) of streamflow"),
             note=TRUE, shades=shades)
  map_annotation()
}
for(d in sort(unique(D$decade))) {
  map_base(); map_sebase()
  shades <- choropleth_cov(spCOV, decade=d, x="se.fit_est_f99", cuts=f99CutsSE, rev=TRUE)
  legend_est(gage=setxt1, title=paste0(d," decade\n",setxt1), note=FALSE, shades=shades, itgage=FALSE)
  map_annotation()
}
dev.off()
#-----------------------------------------------------------------------
pdf("Q99p9fit.pdf", useDingbats=FALSE, width=11, height=10)
for(d in sort(unique(D$decade))) {
  map_base()
  choropleth_decade(D, x="f99.9", cuts=f99p9Cuts, trans=function(t) { log10(t+1) } )
  shades <- choropleth_cov(spCOV, decade=d, x="est_f99.9", cuts=f99p9Cuts)
  legend_est(gage="log(Q99.9+1) of streamflow", title=paste0(d," decade\n","log10(Q99.9+1) of streamflow"),
             note=TRUE, shades=shades)
  map_annotation()
}
for(d in sort(unique(D$decade))) {
  map_base(); map_sebase()
  shades <- choropleth_cov(spCOV, decade=d, x="se.fit_est_f99.9", cuts=f99p9CutsSE, rev=TRUE)
  legend_est(gage=setxt1, title=paste0(d," decade\n",setxt1), note=FALSE, shades=shades, itgage=FALSE)
  map_annotation()
}
dev.off()


plotlmrdia(lmrdia(), xlim=c(.1,1), ylim=c(0,1))
points(spCOV$est_T3[spCOV$decade == "1950"],
       spCOV$est_T4[spCOV$decade == "1950"], col=rgb(1,0,0,.1), cex=0.3, pch=16)
points(spCOV$est_T3[spCOV$decade == "1960"],
       spCOV$est_T4[spCOV$decade == "1960"], col=rgb(0,1,0,.1), cex=0.3, pch=16)
points(spCOV$est_T3[spCOV$decade == "1970"],
       spCOV$est_T4[spCOV$decade == "1970"], col=rgb(0,0,1,.1), cex=0.3, pch=16)
points(spCOV$est_T3[spCOV$decade == "1980"],
       spCOV$est_T4[spCOV$decade == "1980"], col=rgb(1,.5,.5,.1), cex=0.3, pch=17)
points(spCOV$est_T3[spCOV$decade == "1990"],
       spCOV$est_T4[spCOV$decade == "1990"], col=rgb(.5,1,.5,.1), cex=0.3, pch=17)
points(spCOV$est_T3[spCOV$decade == "2000"],
       spCOV$est_T4[spCOV$decade == "2000"], col=rgb(0,.5,1,.1), cex=0.3, pch=17)

write.table(spCOV, file="spCOV.txt", sep=",", quote=TRUE)
save(spCOV, file="spCOV.RData")
write_feather(slot(spCOV, "data"), path="gam_estimates_huc12.feather")

