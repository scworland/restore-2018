library(feather)
library(mgcv)
library(sp)

LATLONG <- paste0("+init=epsg:4269 +proj=longlat +ellps=GRS80 ",
                  "+datum=NAD83 +no_defs +towgs84=0,0,0")
LATLONG <- sp::CRS(LATLONG)
ALBEA <- paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 ",
                "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
ALBEA <- sp::CRS(ALBEA)

east_grids  <- seq(80,100,by=2)
north_grids <- seq(28,38, by=2)

gx <- gy <- vector(mode="numeric")
for(i in 1:length(north_grids)) {
   gy <- c(gy, rep(north_grids[i], length(east_grids)))
   gx <- c(gx, east_grids)
}
GL <- SpatialPoints(cbind(-gx, gy), proj4string=LATLONG)
GL <- spTransform(GL, ALBEA)
XY <- coordinates(GL)
x <- XY[,1]/1000
y <- XY[,2]/1000
ind <- mgcv::inSide(bnd,x,y)
XY <- XY[ind,]
x <- XY[,1]
y <- XY[,2]
plot(x, y, col=2)
GL <- SpatialPoints(cbind(x,y), proj4string=ALBEA)









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



#load(file.choose()) # spDNI_1998to2009.RData
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

message("REMOVING nodata Bed Permeability")
length(spCOV$comid[spCOV$bedperm == "nodata"]) # [1] 198
spCOV <- spCOV[spCOV$bedperm != "nodata",]

EPo <- predict(GM2, newdata=spCOV, se.fit=TRUE)
EP <- EPo$fit; EPo$se.fit[is.na(EPo$se.fit)] <- max(EPo$se.fit, na.rm=TRUE)
EP[EP > log10(3653)] <- log10(3653)
pplo_est <- (3653-10^EP)/3653
spCOV$pplo_est <- pplo_est
spCOV$pplo_est[is.na(spCOV$pplo_est)] <- 0
spCOV$pplo_est[spCOV$pplo_est < 1/3653] <- 0

spCOV$pplo_est <- spCOV$pplo_est <- EP
spCOV$pplo_est[is.na(spCOV$pplo_est)] <- mean(EP, na.rm=TRUE)


quantile(spCOV$pplo_est, probs=(1:9)/10, na.rm=TRUE)
quantile(EPo$se.fit,    probs=(1:9)/10, na.rm=TRUE)

myCuts <- function(x, n=9, ...) {
   labs <- 1:n
   #cuts <- c(0, 0.005, .01, 0.04, .05, 0.1, 0.2, 0.40, 0.60)
   cuts <- c(1.2, 1.4, 1.6, 1.8, 1.9, 2.1, 2.3, 2.6, 3.14)
   cuts <- cuts[labs]
   names(cuts) <- paste("#", labs, sep = "")
   cuts
}

myCuts2 <- function(x, n=9, ...) {
   labs <- 1:n
   #cuts <- c(0.01, 0.015, 0.017, 0.018, 0.019, 0.020, 0.021, 0.023, 0.026)
   cuts <- c(0.015, 0.016, 0.017, 0.019, 0.020, 0.021, 0.022, 0.024, 0.028)
   cuts <- cuts[labs]
   names(cuts) <- paste("#", labs, sep = "")
   cuts
}

library(GISTools)
library(RColorBrewer)



# load(file.choose()) # "GulfStates.RData"
pdf("L1fit.pdf", useDingbats=FALSE, width=11, height=9.5)
  par(lend=1, ljoin=1)
  plot(spCOV, pch=NA); plot(GulfStates, add=TRUE, lty=0, col=grey(0.95))
  polygon(bnd[[1]]$x*1000,bnd[[1]]$y*1000, col=grey(1), lwd=.7)
  tmp <- log10(Z$L1)
  shades <- auto.shading(tmp, cutter=myCuts, n=9,
                         cols=add.alpha(rev(brewer.pal(10,"Spectral")),.7))
  choropleth(Z, tmp, pch=1, lwd=0.7, cex=0.8, shading=shades, add=TRUE)

  tmp <- spCOV$pplo_est; #
  shades <- auto.shading(tmp, cutter=myCuts, n=9,
                         cols=add.alpha(rev(brewer.pal(10,"Spectral")),.2))
  choropleth(spCOV, tmp, pch=16, cex=0.4, shading=shades, add=TRUE)
  shades <- auto.shading(tmp, cutter=myCuts, n=9, cols=rev(brewer.pal(10,"Spectral")))

  ss <- list(x=-100000, y=420000)
  SpatialPolygonsRescale(layout.scale.bar(height=.07), offset=ss, scale=200*1000,
                         fill=c("transparent", "black"), plot.grid=FALSE, lwd=0.4)
  ss <- list(x=-100000, y=403000)
  SpatialPolygonsRescale(layout.scale.bar(height=.08), offset=ss, scale=100*1609.344,
                         fill=c("transparent", "black"), plot.grid=FALSE, lwd=0.4)
  xx <- -142000
  sl <- list(x=xx, y=450000)
  sr <- list(x=xx+200*1000, y=450000)
  text(sl, "0",xx, cex=0.6, pos=4); text(sr, "200 kilometers", cex=0.5, pos=4)
  sl <- list(x=xx, y=383000)
  sr <- list(x=xx+100*1609.344, y=383000)
  text(sl, "0", cex=0.6, pos=4); text(sr, "100 miles", cex=0.5, pos=4)

  txt <- paste0("Albers Equal Area Projection\n",
                "North American Datum of 1983\n",
                "Base modified from U.S. Geological Survey digital data, 1:24,000")
  text(-130000, 330000, txt, pos=4, cex=0.35)
  plot(GulfStates, add=TRUE, lwd=.4, lty=2)
  STATES <- c("Texas", "Oklahoma", "Missouri", "Arkansas", "Louisiana", "Mississippi",
              "Tennessee", "Kentucky", "Alabama", "Georgia", "Florida")
  STATES <- data.frame(easting=c(-440000, -202900.4,  178000,  178000,  279961.5,
                                  430000,  690000,  690000,  690000,
                                 1100000, 1270000),
                       northing=c(955139.0, 1359890.9, 1558716.4, 1400000,  795368.5,
                                  1170000, 1450000, 1600000, 1325000,
                                  1170000, 800000),
                       state=STATES)
  text(STATES$easting, STATES$northing, STATES$state, pos=4, cex=0.7, col=grey(0.3))
  legend(40000, 625000, c(paste0("U.S. Geological Survey streamflow-gaging station:\n",
                                 "Colored* by observed decadal fraction of no flow conditions"),
                          paste0("COMID for HUC12 of National Hydrography Dataset version 2:\n",
                                 "Colored by estimated decadal fraction of no flow conditions")),
         bty="n", cex=0.8, pt.cex=c(0.8,0.6), lwd=c(0.7,1), lty=c(0,0), pch=c(1,16),
         col=c("#3288BDE6","#D53E4FE6"))
  choro.legend(900000, 780000, shades, cex=0.8, bty="n", box.col=grey(1), bg=grey(1),
               fmt="%g", xjust=0, title="Fraction decadal no flow")
  text(205401.9, 490029.4, "* Note that ", cex=.45)
  dev.off()

