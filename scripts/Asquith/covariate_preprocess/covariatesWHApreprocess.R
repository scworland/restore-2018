library(feather)
library(mgcv)
library(sp)

LATLONG <- paste0("+init=epsg:4269 +proj=longlat +ellps=GRS80 ",
                  "+datum=NAD83 +no_defs +towgs84=0,0,0")
LATLONG <- sp::CRS(LATLONG)
ALBEA <- paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 ",
                "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
ALBEA <- sp::CRS(ALBEA)

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

spCOV$developed <- 2*asin(sqrt(spCOV$developed/100))

message("REMOVING nodata Bed Permeability")
length(spCOV$comid[spCOV$bedperm == "nodata"]) # [1] 198
spCOV <- spCOV[spCOV$bedperm != "nodata",]

EP <- EPo <- predict(GM2, newdata=spCOV)
EP[EP > log10(3653)] <- log10(3653)
pplo_est <- (3653-10^EP)/3653
spCOV$pplo_est <- pplo_est
spCOV$pplo_est[is.na(spCOV$pplo_est)] <- 0
spCOV$pplo_est[spCOV$pplo_est < 1/3653] <- 0

quantile(spCOV$pplo_est, probs=(1:9)/10, na.rm=TRUE)

myCuts <- function(x, n=9, ...) {
   labs <- 1:n
   cuts <- c(0, .005, 0.01, .02, .05, 0.1, 0.2, 0.30, 0.40)
   cuts <- cuts[labs]
   names(cuts) <- paste("#", labs, sep = "")
   cuts
}

library(GISTools)
library(RColorBrewer)

# load(file.choose()) # "GulfStates.RData"
pdf("junk.pdf", useDingbats=FALSE)
plot(spCOV, pch=NA)
plot(GulfStates, add=TRUE)
tmp <- Z$pplo
shades <- auto.shading(tmp, cutter=myCuts, n=9, cols=add.alpha(rev(brewer.pal(10,"Spectral")),.5))
choropleth(Z, tmp, pch=1, lwd=0.7, cex=0.6, bg=8, shading=shades, add=TRUE)

  tmp <- spCOV$pplo_est
  shades <- auto.shading(tmp, cutter=myCuts, n=9, cols=add.alpha(rev(brewer.pal(10,"Spectral")),.2))
#  plot(spCOV, lwd=0.3, cex=.3, pch=NA)
  choropleth(spCOV, tmp, pch=16, cex=0.35, shading=shades, add=TRUE)
  shades <- auto.shading(tmp, cutter=myCuts, n=9, cols=rev(brewer.pal(10,"Spectral")),.2)
  choro.legend(600000, 800000, shades, cex=0.4, bty="o", box.col=grey(1), bg=grey(1))
dev.off()
