library(feather)
library(sp)

load("../../../gis/solar/spDNI_1998to2009.RData") # spDNI_1998to2009.RData
FDC <- read_feather("../../../data/gage/all_gage_covariates.feather") # "all_gage_covariates.feather"
COV <- read_feather("../../../data/huc12/all_huc12_covariates.feather") # "all_huc12_covariates.feather"

sites <- unique(FDC$site_no)
sitefile <- dataRetrieval::readNWISsite(sites)
sitefile <- sitefile[sitefile$agency_cd != "USCE",]
FDC$dec_lat_va <-  NA
FDC$dec_long_va <- NA
for(site in FDC$site_no) {
  FDC$dec_lat_va[FDC$site_no == site] <- sitefile$dec_lat_va[sitefile$site_no == site]
  FDC$dec_long_va[FDC$site_no == site] <- sitefile$dec_long_va[sitefile$site_no == site]
}
summary(abs(FDC$dec_lat_va  - FDC$dec_lat_va ))
summary(abs(FDC$dec_long_va - FDC$dec_long_va))

LATLONG <- paste0("+proj=longlat +ellps=GRS80 ",
                  "+datum=NAD83 +no_defs +towgs84=0,0,0")
LATLONG <- sp::CRS(LATLONG)
ALBEA <- paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 ",
                "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
ALBEA <- sp::CRS(ALBEA)

nFDC <- data.frame(comid=FDC$comid,
                   site_no=FDC$site_no,
                   huc12=FDC$huc12,
                   dec_long_va=FDC$dec_long_va,
                   dec_lat_va=FDC$dec_lat_va,
                   decade=FDC$decade,
                   stringsAsFactors=FALSE)
nFDC <- SpatialPointsDataFrame(cbind(nFDC$dec_long_va, nFDC$dec_lat_va), nFDC,
                               proj4string=LATLONG)
nFDC <- spTransform(nFDC, ALBEA)
SO <- over(nFDC, spDNI_1998to2009)
nFDC$dni_ann <- SO$ANN_DNI
nFDC$dni_jan <- SO$JAN; nFDC$dni_feb <- SO$FEB
nFDC$dni_mar <- SO$MAR; nFDC$dni_apr <- SO$APR
nFDC$dni_may <- SO$MAY; nFDC$dni_jun <- SO$JUN
nFDC$dni_jul <- SO$JUL; nFDC$dni_aug <- SO$AUG
nFDC$dni_sep <- SO$SEP; nFDC$dni_oct <- SO$OCT
nFDC$dni_nov <- SO$NOV; nFDC$dni_dec <- SO$DEC
rm(SO)

nCOV <- data.frame(comid=COV$comid,
                   huc12=COV$huc12,
                   dec_long_va=COV$dec_long_va,
                   dec_lat_va=COV$dec_lat_va,
                   decade=COV$decade,
                   stringsAsFactors=FALSE)
nCOV <- SpatialPointsDataFrame(cbind(nCOV$dec_long_va,nCOV$dec_lat_va), data=nCOV,
                                proj4string=LATLONG)
nCOV <- spTransform(nCOV, ALBEA)
SO <- over(nCOV, spDNI_1998to2009)
nCOV$dni_ann <- SO$ANN_DNI
nCOV$dni_jan <- SO$JAN; nCOV$dni_feb <- SO$FEB
nCOV$dni_mar <- SO$MAR; nCOV$dni_apr <- SO$APR
nCOV$dni_may <- SO$MAY; nCOV$dni_jun <- SO$JUN
nCOV$dni_jul <- SO$JUL; nCOV$dni_aug <- SO$AUG
nCOV$dni_sep <- SO$SEP; nCOV$dni_oct <- SO$OCT
nCOV$dni_nov <- SO$NOV; nCOV$dni_dec <- SO$DEC
rm(SO)

write_feather(slot(nFDC, "data"), "all_gage_solar.feather" )
write_feather(slot(nCOV, "data"), "all_huc12_solar.feather")
