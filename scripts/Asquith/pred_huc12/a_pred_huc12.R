library(feather)
library(mgcv)
library(sp)
library(GISTools)
library(RColorBrewer)
library(lmomco)
source("../fdc_model/gamIntervals.R")

load("../fdc_model/Models.RData")
load("../../../../GIS/GulfStates.RData") # GulfStates.RData
load("../../../../GIS/RESTORE_MGCV_BND.RData") # "RESTORE_MGCV_BND.RData"

LATLONG <- paste0("+proj=longlat +ellps=GRS80 ",
                  "+datum=NAD83 +no_defs +towgs84=0,0,0")
LATLONG <- sp::CRS(LATLONG)
ALBEA <- paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 ",
                "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
ALBEA <- sp::CRS(ALBEA)

COV <- COVo <- read_feather("../../../data/huc12/all_huc12_covariates.feather")
length(COVo$comid)                           # [1] 55128
SO <- read_feather("../../../data/huc12/all_huc12_solar.feather")
COV <- cbind(COV, SO)

spCOV <- SpatialPointsDataFrame(cbind(COV$dec_long_va,COV$dec_lat_va), data=COV,
                                proj4string=LATLONG)
spCOV <- spTransform(spCOV, ALBEA)
XY <- coordinates(spCOV)
spCOV$x <- spCOV$east <- XY[,1]/1000; spCOV$y <- spCOV$north <- XY[,2]/1000
rm(COV, XY)

length(spCOV$nid_storage[ spCOV$basin_area == 0]) # [1] 0
length(spCOV$norm_storage[spCOV$basin_area == 0]) # [1] 0


length(spCOV$flood_storage[spCOV$flood_storage > 1]) # [1] 35
#boxplot(spCOV$flood_storage)
spCOV$flood_storage[spCOV$flood_storage > 1] <- 1 # truncate predictions in the extrapolation zone
#spCOV <- spCOV[spCOV$flood_storage <= 1,] # This is about 1 meter of watershed depth equivalent
# which is close to the maximum observed in the streamgage network itself. But if the GAM
# models PPLO-->T6 don't use flood_storage, just don't delete those.

length(spCOV$comid[spCOV$basin_area == 0]) # [1] 0

spCOV$ppt_mean       <- log10(spCOV$ppt_mean)
spCOV$temp_mean      <- log10(spCOV$temp_mean)
spCOV$basin_area     <- log10(spCOV$basin_area)
spCOV$basin_slope    <- log10(spCOV$basin_slope/100)
# DD is the "data" and not the prediction points.
flood_storage_offset <- 1E-6; # first even log10 cycle below min(DD$flood_storage[DD$flood_storage > 0])
spCOV$flood_storage <- log10(spCOV$flood_storage + flood_storage_offset)

length(spCOV$comid[! is.finite(spCOV$ppt_mean)])     # [1] 0
length(spCOV$comid[! is.finite(spCOV$temp_mean)])    # [1] 0
length(spCOV$comid[! is.finite(spCOV$basin_area)])   # [1] 0
length(spCOV$comid[! is.finite(spCOV$basin_slope)])  # [1] 0

spCOV$comid.1       <- NULL
spCOV$huc12.1       <- NULL
spCOV$dec_long_va.1 <- NULL
spCOV$dec_lat_va.1  <- NULL
spCOV$decade.1      <- NULL

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
unique(spCOV$site_no[spCOV$ed_rch_zone == "1"]) # should be zero as we decide not to reintersect?

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


quantile(spCOV$grassland,    probs=(1:9)/10, na.rm=TRUE)

grassCuts <- function(x, n=9, ...) {
   labs <- 1:n
   cuts <- c(1, 5, 10, 15, 25, 30, 35, 40, 45)/100
   cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}

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
   cuts <- c(0.0162, 0.0177, 0.0200, 0.0228, 0.0254, 0.0276, 0.0302, 0.0340, 0.0422)
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
  cuts <- c(-0.4, -0.1, 0.08, 0.22, 0.37, 0.55, 0.76, 1.07, 1.58)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}

L1CutsSE <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.0157, 0.0165, 0.0172, 0.0178, 0.0185, 0.0193, 0.0203, 0.02149, 0.0250)
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
H12T2df$rse_T2 <- ET2$residual.scale[1]
H12T2df$se.fit_T2 <- ET2$se.fit
H12T2df <- SpatialPointsDataFrame(cbind(H12T2df$dec_long_va,H12T2df$dec_lat_va),
                                    data=H12T2df, proj4string=LATLONG)
H12T2df <- spTransform(H12T2df, ALBEA)


quantile(H12T2df$est_T2,    probs=(1:9)/10, na.rm=TRUE)
quantile(H12T2df$se.fit_T2, probs=(1:9)/10, na.rm=TRUE)

T2Cuts <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.54, 0.60, 0.63, 0.67, 0.70, 0.73, 0.77, 0.81, 0.86)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}

T2CutsSE <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.0106, 0.0111, 0.0115, 0.0120, 0.0125, 0.0130, 0.0138, 0.0150, 0.0170)
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
H12T3df$rse_T3 <- ET3$residual.scale[1]
H12T3df$se.fit_T3 <- ET3$se.fit
H12T3df <- SpatialPointsDataFrame(cbind(H12T3df$dec_long_va,H12T3df$dec_lat_va),
                                    data=H12T3df, proj4string=LATLONG)
H12T3df <- spTransform(H12T3df, ALBEA)


quantile(H12T3df$est_T3,    probs=(1:9)/10, na.rm=TRUE)
quantile(H12T3df$se.fit_T3, probs=(1:9)/10, na.rm=TRUE)

T3Cuts <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.53, 0.57, 0.60, 0.63, 0.66, 0.68, 0.71, 0.77, 0.84)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}

T3CutsSE <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.0100, 0.0103, 0.0107, 0.0110, 0.0115, 0.0120, 0.0127, 0.0138, 0.0158)
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
H12T4df$rse_T4 <- ET4$residual.scale[1]
H12T4df$se.fit_T4 <- ET4$se.fit
H12T4df <- SpatialPointsDataFrame(cbind(H12T4df$dec_long_va,H12T4df$dec_lat_va),
                                    data=H12T4df, proj4string=LATLONG)
H12T4df <- spTransform(H12T4df, ALBEA)

quantile(H12T4df$est_T4,    probs=(1:9)/10, na.rm=TRUE)
quantile(H12T4df$se.fit_T4, probs=(1:9)/10, na.rm=TRUE)

T4Cuts <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.33, 0.39, 0.42, 0.45, 0.47, 0.50, 0.53, 0.62, 0.72)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}

T4CutsSE <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.0105, 0.0110, 0.0115, 0.0119, 0.0123, 0.0129, 0.0137, 0.0149, 0.0171)
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
H12T5df$rse_T5 <- ET5$residual.scale[1]
H12T5df$se.fit_T5 <- ET5$se.fit
H12T5df <- SpatialPointsDataFrame(cbind(H12T5df$dec_long_va,H12T5df$dec_lat_va),
                                    data=H12T5df, proj4string=LATLONG)
H12T5df <- spTransform(H12T5df, ALBEA)

quantile(H12T5df$est_T5,    probs=(1:9)/10, na.rm=TRUE)
quantile(H12T5df$se.fit_T5, probs=(1:9)/10, na.rm=TRUE)


T5Cuts <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.22, 0.28, 0.32, 0.34, 0.36, 0.39, 0.42, 0.50, 0.61)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}

T5CutsSE <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(0.0102, 0.0107, 0.0111, 0.0115, 0.0119,
            0.0125, 0.0132, 0.0144, 0.0164)
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

# source("a_plot_huc12.R")
