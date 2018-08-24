library(feather)
library(mgcv)
library(sp)
library(GISTools)
library(RColorBrewer)
library(lmomco)
source("a_basemap_funcs.R")
source("../fdc_model/gamIntervals.R")

load("../fdc_model/FDCEST.RData")
load("../fdc_model/Models.RData")
load("../../../../GIS/GulfStates.RData") # GulfStates.RData
#load(file.choose()) # spDNI_1998to2009.RData
load("../../../../GIS/RESTORE_MGCV_BND.RData") # "RESTORE_MGCV_BND.RData"

LATLONG <- paste0("+proj=longlat +ellps=GRS80 ",
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
x <- XY[,1]; y <- XY[,2]
#ind <- mgcv::inSide(bnd,x,y)
#XY <- XY[ind,]
GL <- SpatialPoints(cbind(x,y), proj4string=ALBEA)
ix <- 1:length(x)
plot(GL, pch=1, col=2)
text(XY[,1],XY[,2], ix)
GL <- GL[-c(1,3:9, 12, 14:20, 23, 34, 45)]

COV <- COVo <- read_feather("../../../data/huc12/all_huc12_covariates.feather")
length(COVo$comid)                           # [1] 55128
SO <- read_feather("../../../data/huc12/all_huc12_solar.feather")
COV <- cbind(COV, SO)

spCOV <- SpatialPointsDataFrame(cbind(COV$lon,COV$lat), data=COV,
                                proj4string=LATLONG)
spCOV <- spTransform(spCOV, ALBEA)
XY <- coordinates(spCOV)
spCOV$x <- spCOV$east <- XY[,1]/1000; spCOV$y <- spCOV$north <- XY[,2]/1000

rm(COV, XY)
#SO <- over(spCOV, spDNI_1998to2009)

length(spCOV$nid_storage[ spCOV$basin_area == 0])
length(spCOV$norm_storage[spCOV$basin_area == 0])


length(spCOV$flood_storage[spCOV$flood_storage > 1]) # [1] 35
boxplot(spCOV$flood_storage)
#spCOV <- spCOV[spCOV$flood_storage <= 1,] # This is about 1 meter of watershed depth equivalent
# which is close to the maximum observed in the streamgage network itself. But if the GAM
# models PPLO-->T6 don't use flood_storage, just don't delete those.

length(spCOV$comid[spCOV$basin_area == 0])

spCOV$ppt_mean        <- log10(spCOV$ppt_mean)
spCOV$temp_mean       <- log10(spCOV$temp_mean)
spCOV$basin_area  <- log10(spCOV$basin_area)
spCOV$basin_slope <- log10(spCOV$basin_slope/100)
length(spCOV$comid[! is.finite(spCOV$ppt_mean)])     # [1] 0
length(spCOV$comid[! is.finite(spCOV$temp_mean)])    # [1] 0
length(spCOV$comid[! is.finite(spCOV$basin_area)])   # [1] 0
length(spCOV$comid[! is.finite(spCOV$basin_slope)])  # [1] 0


spCOV$decade   <- as.factor(spCOV$decade)
spCOV$bedperm  <- as.factor(spCOV$bedperm)
spCOV$aquifers <- as.factor(spCOV$aquifers)
spCOV$soller   <- as.factor(spCOV$soller)
spCOV$hlr      <- as.factor(spCOV$hlr)
spCOV$ecol3    <- as.factor(spCOV$ecol3)
spCOV$physio   <- as.factor(spCOV$physio)
spCOV$statsgo  <- as.factor(spCOV$statsgo)


spCOV$decade       <- as.factor(spCOV$decade);       levels(spCOV$decade)
spCOV$cat_soller   <- as.factor(spCOV$cat_soller);   levels(spCOV$cat_soller)
spCOV$soller       <- as.factor(spCOV$soller);       levels(spCOV$soller)
spCOV$cat_aquifers <- as.factor(spCOV$cat_aquifers); levels(spCOV$cat_aquifers)
spCOV$aquifers     <- as.factor(spCOV$aquifers);     levels(spCOV$aquifers)
spCOV$bedperm      <- as.factor(spCOV$bedperm);      levels(spCOV$bedperm)
spCOV$cat_physio   <- as.factor(spCOV$cat_physio);   levels(spCOV$cat_physio)
spCOV$physio       <- as.factor(spCOV$physio);       levels(spCOV$physio)
spCOV$cat_ecol3    <- as.factor(spCOV$cat_ecol3);    levels(spCOV$cat_ecol3)
spCOV$ecol3        <- as.factor(spCOV$ecol3);        levels(spCOV$ecol3)
spCOV$hlr          <- as.factor(spCOV$hlr);          levels(spCOV$hlr)
spCOV$statsgo      <- as.factor(spCOV$statsgo);      levels(spCOV$statsgo)
spCOV$ed_rch_zone  <- as.factor(spCOV$ed_rch_zone);  levels(spCOV$ed_rch_zone)
unique(spCOV$site_no[spCOV$ed_rch_zone == "1"])

summary(spCOV$cat_soller)
summary(spCOV$soller)
summary(spCOV$cat_aquifers)
summary(spCOV$aquifers)
summary(spCOV$bedperm)
summary(spCOV$cat_physio)
summary(spCOV$physio)
summary(spCOV$cat_ecol3)
summary(spCOV$ecol3)
summary(spCOV$hlr)
summary(spCOV$statsgo)
summary(spCOV$ed_rch_zone)


# Transformation and Retransformation Functions for the Sin Transformation of
# Percentile data
dotransin <- function(p) 2*asin(sqrt(p/100))
retransin <- function(p)    sin(p/2)^2*100

spCOV$barren              <-  dotransin(spCOV$barren)
spCOV$cultivated_cropland <-  dotransin(spCOV$cultivated_cropland)
spCOV$deciduous_forest    <-  dotransin(spCOV$deciduous_forest)
spCOV$developed           <-  dotransin(spCOV$developed)
spCOV$evergreen_forest    <-  dotransin(spCOV$evergreen_forest)
spCOV$grassland           <-  dotransin(spCOV$grassland)
spCOV$hay_pasture         <-  dotransin(spCOV$hay_pasture)
spCOV$herbaceous_wetland  <-  dotransin(spCOV$herbaceous_wetland)
spCOV$mixed_forest        <-  dotransin(spCOV$mixed_forest)
spCOV$shrubland           <-  dotransin(spCOV$shrubland)
spCOV$water               <-  dotransin(spCOV$water)
spCOV$woody_wetland       <-  dotransin(spCOV$woody_wetland)
spCOV$bfi                 <-  dotransin(spCOV$bfi)



EPo <- predict(PPLO, newdata=spCOV, se.fit=TRUE); length(EPo$fit)
EPo <- gamIntervals(EPo, gam=PPLO, interval="prediction", sigma=PPLO$pplo.sigma[1])

# Terms invert in upper/lower meaning and hence the flipping during data.frame construction.
H12PPLOdf <- data.frame(comid=spCOV$comid, huc12=spCOV$huc12,
                     decade=as.character(spCOV$decade),
                     dec_long_va=spCOV$dec_long_va, dec_lat_va=spCOV$dec_lat_va,
                     est_lwr_pplo=(3653-10^EPo$upr)/3653,
                     est_pplo    =(3653-10^EPo$fit)/3653,
                     est_upr_pplo=(3653-10^EPo$lwr)/3653,
                     est_lwr_flowtime=EPo$lwr,
                     est_flowtime=EPo$fit,
                     est_upr_flowtime=EPo$upr,
                     stringsAsFactors=FALSE)
H12PPLOdf$est_lwr_pplo[H12PPLOdf$est_lwr_pplo < 0] <- 0
H12PPLOdf$est_pplo[    H12PPLOdf$est_pplo     < 0] <- 0
H12PPLOdf$est_upr_pplo[H12PPLOdf$est_upr_pplo < 0] <- 0
H12PPLOdf$rse_pplo <- PPLO$pplo.sigma[1]
H12PPLOdf$se.fit_pplo <- EPo$se.fit
H12PPLOdf <- SpatialPointsDataFrame(cbind(H12PPLOdf$dec_long_va,H12PPLOdf$dec_lat_va),
                                    data=H12PPLOdf, proj4string=LATLONG)
H12PPLOdf <- spTransform(H12PPLOdf, ALBEA)


quantile(H12PPLOdf$est_pplo,    probs=(1:9)/10, na.rm=TRUE)
quantile(H12PPLOdf$se.fit_pplo, probs=(1:9)/10, na.rm=TRUE)

pploCuts <- function(x, n=9, ...) {
   labs <- 1:n
   cuts <- c(0, 0.005, .01, 0.04, .05, 0.1, 0.2, 0.4, 0.6)
   cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}

pploCutsSE <- function(x, n=9, ...) {
   labs <- 1:n
   cuts <- c(0.0152, 0.0167, 0.0189, 0.0213, 0.0236, 0.0253, 0.0274, 0.0301, 0.0354)
   cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}

EL1 <- predict(L1, newdata=spCOV, se.fit=TRUE)
EL1 <- gamIntervals(EL1, gam=L1, interval="prediction")
H12L1df <- data.frame(comid=spCOV$comid, huc12=spCOV$huc12,
                     decade=spCOV$decade,
                     dec_long_va=spCOV$dec_long_va, dec_lat_va=spCOV$dec_lat_va,
                     bias_corr=L1$duan_smearing,
                     est_lwr_L1=10^EL1$lwr,
                     est_L1    =10^EL1$fit,
                     est_upr_L1=10^EL1$upr, stringsAsFactors=FALSE)
H12L1df$rse_L1 <- EL1$residual.scale[1]
H12L1df$se.fit_L1 <- EL1$se.fit
H12L1df <- SpatialPointsDataFrame(cbind(H12L1df$dec_long_va,H12L1df$dec_lat_va),
                                    data=H12L1df, proj4string=LATLONG)
H12L1df <- spTransform(H12L1df, ALBEA)


quantile(log10(H12L1df$est_L1), probs=(1:9)/10, na.rm=TRUE)
quantile(H12L1df$se.fit_L1,     probs=(1:9)/10, na.rm=TRUE)

L1Cuts <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(-0.383, -0.104, 0.067, 0.219, 0.375, 0.546, 0.758, 1.071, 1.575)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}

L1CutsSE <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.0157, 0.0165, 0.0171, 0.0177, 0.0183, 0.019, 0.020, 0.0214, 0.0243)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}



ET2 <- predict(T2, newdata=spCOV, se.fit=TRUE)
ET2 <- gamIntervals(ET2, gam=T2, interval="prediction")
H12T2df <- data.frame(comid=spCOV$comid, huc12=spCOV$huc12,
                     decade=spCOV$decade,
                     dec_long_va=spCOV$dec_long_va, dec_lat_va=spCOV$dec_lat_va,
                     est_lwr_T2=ET2$lwr,
                     est_T2    =ET2$fit,
                     est_upr_T2=ET2$upr, stringsAsFactors=FALSE)
H12T2df$rse_T2 <- sigma
H12T2df$se.fit_T2 <- ET2$se.fit
H12T2df <- SpatialPointsDataFrame(cbind(H12T2df$dec_long_va,H12T2df$dec_lat_va),
                                    data=H12T2df, proj4string=LATLONG)
H12T2df <- spTransform(H12T2df, ALBEA)


quantile(H12T2df$est_T2,    probs=(1:9)/10, na.rm=TRUE)
quantile(H12T2df$se.fit_T2, probs=(1:9)/10, na.rm=TRUE)

T2Cuts <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.55, 0.60, 0.64, 0.67, 0.70, 0.73, 0.76, 0.80, 0.86)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}

T2CutsSE <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.0106, 0.0111, 0.0115, 0.0119, 0.0123, 0.0127, 0.0134, 0.0144, 0.0163)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}


ET3 <- predict(T3, newdata=spCOV, se.fit=TRUE)
ET3 <- gamIntervals(ET3, gam=T3, interval="prediction")
H12T3df <- data.frame(comid=spCOV$comid, huc12=spCOV$huc12,
                     decade=spCOV$decade,
                     dec_long_va=spCOV$dec_long_va, dec_lat_va=spCOV$dec_lat_va,
                     est_lwr_T3=ET3$lwr,
                     est_T3    =ET3$fit,
                     est_upr_T3=ET3$upr, stringsAsFactors=FALSE)
H12T3df$rse_T3 <- sigma
H12T3df$se.fit_T3 <- ET3$se.fit
H12T3df <- SpatialPointsDataFrame(cbind(H12T3df$dec_long_va,H12T3df$dec_lat_va),
                                    data=H12T3df, proj4string=LATLONG)
H12T3df <- spTransform(H12T3df, ALBEA)


quantile(H12T3df$est_T3,    probs=(1:9)/10, na.rm=TRUE)
quantile(H12T3df$se.fit_T3, probs=(1:9)/10, na.rm=TRUE)

T3Cuts <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.53, 0.57, 0.60, 0.63, 0.66, 0.68, 0.71, 0.76, 0.83)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}

T3CutsSE <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.0104, 0.0109, 0.0113, 0.0117, 0.0121, 0.0126, 0.0132, 0.0142, 0.0161)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}


ET4 <- predict(T4, newdata=spCOV, se.fit=TRUE)
ET4 <- gamIntervals(ET4, gam=T4, interval="prediction")
H12T4df <- data.frame(comid=spCOV$comid, huc12=spCOV$huc12,
                     decade=spCOV$decade,
                     dec_long_va=spCOV$dec_long_va, dec_lat_va=spCOV$dec_lat_va,
                     est_lwr_T4=ET4$lwr,
                     est_T4    =ET4$fit,
                     est_upr_T4=ET4$upr, stringsAsFactors=FALSE)
H12T4df$rse_T4 <- sigma
H12T4df$se.fit_T4 <- ET4$se.fit
H12T4df <- SpatialPointsDataFrame(cbind(H12T4df$dec_long_va,H12T4df$dec_lat_va),
                                    data=H12T4df, proj4string=LATLONG)
H12T4df <- spTransform(H12T4df, ALBEA)

quantile(H12T4df$est_T4,    probs=(1:9)/10, na.rm=TRUE)
quantile(H12T4df$se.fit_T4, probs=(1:9)/10, na.rm=TRUE)

T4Cuts <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.34, 0.40, 0.43, 0.46, 0.47, 0.50, 0.54, 0.60, 0.71)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}

T4CutsSE <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.0113, 0.0118, 0.0122, 0.0125, 0.0130, 0.0135, 0.0142, 0.0154, 0.0175)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}


ET5 <- predict(T5, newdata=spCOV, se.fit=TRUE)
ET5 <- gamIntervals(ET5, gam=T5, interval="prediction")
H12T5df <- data.frame(comid=spCOV$comid, huc12=spCOV$huc12,
                     decade=spCOV$decade,
                     dec_long_va=spCOV$dec_long_va, dec_lat_va=spCOV$dec_lat_va,
                     est_lwr_T5=ET5$lwr,
                     est_T5    =ET5$fit,
                     est_upr_T5=ET5$upr, stringsAsFactors=FALSE)
H12T5df$rse_T5 <- sigma
H12T5df$se.fit_T5 <- ET5$se.fit
H12T5df <- SpatialPointsDataFrame(cbind(H12T5df$dec_long_va,H12T5df$dec_lat_va),
                                    data=H12T5df, proj4string=LATLONG)
H12T5df <- spTransform(H12T5df, ALBEA)

quantile(H12T5df$est_T5,    probs=(1:9)/10, na.rm=TRUE)
quantile(H12T5df$se.fit_T5, probs=(1:9)/10, na.rm=TRUE)


T5Cuts <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.23, 0.29, 0.32, 0.35, 0.37, 0.39, 0.42, 0.48, 0.60)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}

T5CutsSE <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.0111, 0.0116, 0.0120, 0.0124, 0.0128,
            0.0133, 0.0140, 0.0151, 0.0172)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}


#-----------------------------------------------------------------------

H12PPLOdf$delta_est_pplo <- NA
for(comid in unique(H12PPLOdf$comid)) {
  H12PPLOdf$delta_est_pplo[        H12PPLOdf$comid == comid] <-
      c(NA,diff(H12PPLOdf$est_pplo[H12PPLOdf$comid == comid]))
}
H12L1df$delta_est_L1 <- NA
for(comid in unique(H12L1df$comid)) {
   H12L1df$delta_est_L1[       H12L1df$comid == comid] <-
      c(NA,diff(H12L1df$est_L1[H12L1df$comid == comid]))
}
H12T2df$delta_est_T2 <- NA
for(comid in unique(H12T2df$comid)) {
   H12T2df$delta_est_T2[       H12T2df$comid == comid] <-
      c(NA,diff(H12T2df$est_T2[H12T2df$comid == comid]))
}
H12T3df$delta_est_T3 <- NA
for(comid in unique(H12T3df$comid)) {
   H12T3df$delta_est_T3[       H12T3df$comid == comid] <-
      c(NA,diff(H12T3df$est_T3[H12T3df$comid == comid]))
}
H12T4df$delta_est_T4 <- NA
for(comid in unique(H12T4df$comid)) {
   H12T4df$delta_est_T4[       H12T4df$comid == comid] <-
      c(NA,diff(H12T4df$est_T4[H12T4df$comid == comid]))
}
H12T5df$delta_est_T5 <- NA
for(comid in unique(H12T5df$comid)) {
   H12T5df$delta_est_T5[       H12T5df$comid == comid] <-
      c(NA,diff(H12T5df$est_T5[H12T5df$comid == comid]))
}

write_feather(slot(H12PPLOdf, "data"), "all_gam_huc12_pplo.feather")
write_feather(slot(H12L1df,   "data"),   "all_gam_huc12_L1.feather")
write_feather(slot(H12T2df,   "data"),   "all_gam_huc12_T2.feather")
write_feather(slot(H12T3df,   "data"),   "all_gam_huc12_T3.feather")
write_feather(slot(H12T4df,   "data"),   "all_gam_huc12_T4.feather")
write_feather(slot(H12T5df,   "data"),   "all_gam_huc12_T5.feather")



#-----------------------------------------------------------------------
pdf("PPLOfit_junk.pdf", useDingbats=FALSE, width=11, height=10)
  plot(spRESTORE_MGCV_BND)  # by creation of the PDF, we can get a handle on a global
  usr <- par()$usr # setting of the plotting limits by preserving the usr.
dev.off()
unlink("PPLOfit_junk.pdf")  # just quietly throw the file away

pdf("PPLOfit.pdf", useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    map_base(xlim=usr[1:2], ylim=usr[3:4])
    choropleth_decade(D, x="pplo", cuts=pploCuts, rev=TRUE)
    shades <- choropleth_cov(H12PPLOdf, decade=d, x="est_pplo", cuts=pploCuts, rev=TRUE)
    legend_est(gage="no flow fraction", title=paste0(d," decade\n","no flow fraction"),
               note=TRUE, shades=shades)
    map_annotation()
  }
dev.off()
pdf("PPLOsefit.pdf", useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    map_base(xlim=usr[1:2], ylim=usr[3:4]); map_sebase()
    shades <- choropleth_cov(H12PPLOdf, decade=d, x="se.fit_pplo", cuts=pploCutsSE, rev=TRUE)
    legend_est(gage=setxt1, title=paste0(d,"decade\n",setxt1), note=FALSE, shades=shades, itgage=FALSE)
    map_annotation()
  }
dev.off()
#-----------------------------------------------------------------------
pdf("L1fit.pdf", useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    map_base(xlim=usr[1:2], ylim=usr[3:4])
    choropleth_decade(D, x="L1", cuts=L1Cuts, trans=log10)
    shades <- choropleth_cov(H12L1df, decade=d, x="est_L1", cuts=L1Cuts)
    legend_est(gage="mean streamflow", title=paste0(d," decade\n","mean streamflow, in log10(cms)"),
               note=TRUE, shades=shades)
    map_annotation()
  }
dev.off()
pdf("L1sefit.pdf", useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    map_base(xlim=usr[1:2], ylim=usr[3:4]); map_sebase()
    shades <- choropleth_cov(H12L1df, decade=d, x="se.fit_L1", cuts=L1CutsSE, rev=TRUE)
    legend_est(gage=setxt1, title=paste0(d," decade\n",setxt1), note=FALSE, shades=shades, itgage=FALSE)
    map_annotation()
  }
dev.off()
#-----------------------------------------------------------------------
pdf("T2fit.pdf", useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    map_base(xlim=usr[1:2], ylim=usr[3:4]);
    choropleth_decade(D, x="T2", cuts=T2Cuts)
    shades <- choropleth_cov(H12T2df, decade=d, x="est_T2", cuts=T2Cuts)
    legend_est(gage="L-CV of streamflow", title=paste0(d," decade\n","L-CV of streamflow"),
               note=TRUE, shades=shades)
   map_annotation()
  }
dev.off()
pdf("T2sefit.pdf", useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    map_base(xlim=usr[1:2], ylim=usr[3:4]); map_sebase()
    shades <- choropleth_cov(H12T2df, decade=d, x="se.fit_T2", cuts=T2CutsSE, rev=TRUE)
    legend_est(gage=setxt1, title=paste0(d," decade\n",setxt1), note=FALSE, shades=shades, itgage=FALSE)
    map_annotation()
  }
dev.off()
#-----------------------------------------------------------------------
pdf("T3fit.pdf", useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    map_base(xlim=usr[1:2], ylim=usr[3:4])
    choropleth_decade(D, x="T3", cuts=T3Cuts)
    shades <- choropleth_cov(H12T3df, decade=d, x="est_T3", cuts=T3Cuts)
    legend_est(gage="L-skew of streamflow", title=paste0(d," decade\n","L-skew of streamflow"),
               note=TRUE, shades=shades)
    map_annotation()
  }
dev.off()
pdf("T3sefit.pdf", useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    map_base(xlim=usr[1:2], ylim=usr[3:4]); map_sebase()
    shades <- choropleth_cov(H12T3df, decade=d, x="se.fit_T3", cuts=T3CutsSE, rev=TRUE)
    legend_est(gage=setxt1, title=paste0(d," decade\n",setxt1), note=FALSE, shades=shades, itgage=FALSE)
    map_annotation()
  }
dev.off()
#-----------------------------------------------------------------------
pdf("T4fit.pdf", useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    map_base(xlim=usr[1:2], ylim=usr[3:4])
    choropleth_decade(D, x="T4", cuts=T4Cuts)
    shades <- choropleth_cov(H12T4df, decade=d, x="est_T4", cuts=T4Cuts)
    legend_est(gage="L-kurtosis of streamflow", title=paste0(d," decade\n","L-kurtosis of streamflow"),
               note=TRUE, shades=shades)
    map_annotation()
  }
dev.off()
pdf("T4sefit.pdf", useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    map_base(xlim=usr[1:2], ylim=usr[3:4]); map_sebase()
    shades <- choropleth_cov(H12T4df, decade=d, x="se.fit_T4", cuts=T4CutsSE, rev=TRUE)
    legend_est(gage=setxt1, title=paste0(d," decade\n",setxt1), note=FALSE, shades=shades, itgage=FALSE)
    map_annotation()
  }
dev.off()
#-----------------------------------------------------------------------
pdf("T5fit.pdf", useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    map_base(xlim=usr[1:2], ylim=usr[3:4])
    choropleth_decade(D, x="T5", cuts=T5Cuts)
    shades <- choropleth_cov(H12T5df, decade=d, x="est_T5", cuts=T5Cuts)
    legend_est(gage="Tau5 of streamflow", title=paste0(d," decade\n","Tau5 of streamflow"),
               note=TRUE, shades=shades)
    map_annotation()
  }
dev.off()
pdf("T5sefit.pdf", useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    map_base(xlim=usr[1:2], ylim=usr[3:4]); map_sebase()
    shades <- choropleth_cov(H12T5df, decade=d, x="se.fit_T5", cuts=T5CutsSE, rev=TRUE)
    legend_est(gage=setxt1, title=paste0(d," decade\n",setxt1), note=FALSE, shades=shades, itgage=FALSE)
    map_annotation()
  }
dev.off()


quantile(H12L1df$delta_est_L1, probs=(1:9)/10, na.rm=TRUE)

L1delCuts <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(-1, -.5, -.2, -0.05, 0, 0.05, 0.2, 0.5, 1)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}

pdf("L1del.pdf", useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    if(d == "1950") next
    map_base(xlim=usr[1:2], ylim=usr[3:4])
    choropleth_decade(D, x="L1", cuts=L1delCuts)
    shades <- choropleth_cov(H12L1df, decade=d, x="delta_est_L1", cuts=L1delCuts)
    legend_est(gage="L1 of streamflow", title=paste0(d," decade\n","Change in L1 of streamflow"),
               note=TRUE, shades=shades)
    map_annotation()
  }
dev.off()
