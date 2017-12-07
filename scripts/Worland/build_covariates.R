
library(pacman)
pacman::p_load(tidyverse,stringr,lubridate,feather,geoknife,sf,sp,dataRetrieval)
source('scripts/worland/nldi_funs.R')

# NHD+ basin characteristics ----
# load HUC12s for restore footprint
rest_huc12 <- st_read('data/shapefiles/restore_hucs/restoreHUC12s.shp', stringsAsFactors = F)

# subset mobile basin
mobile <- filter(rest_huc12, str_detect(HUC_8, '0316')) %>%
  distinct(HUC_12, .keep_all=T)

st_write(mobile,'data/shapefiles/mobile_shpfiles/mobileHUC12s.shp')

# load gages site list
sites <- read_csv("data/full_site_list.csv") 

# find coordinates of gages in mobile
site_info <- readNWISdata(sites=sites$siteno, service="site") %>%
  select(siteno=site_no,lon=dec_long_va,lat=dec_lat_va) %>%
  st_as_sf(coords=c("lon","lat"), crs=st_crs(rest_huc12), remove=F) %>%
  #st_intersection(mobile) %>%
  select(siteno,lon,lat,HUC_12) %>%
  na.omit() 

st_write(site_info, 'data/shapefiles/mobile_shpfiles/mobile_gages.shp')
st_write(site_info, 'data/shapefiles/gages/restore_gages.shp')

# extract basin characteristics for nwis sites
gage_chars <- NULL
for(i in 1:nrow(site_info)){
  print(i)
  site <- site_info$siteno[i]
  char_list <- characteristics_nldi("ss_gages",site,tier = "test")
  
  if(!is.null(char_list)) {
    temp_df <- char_list$characteristics[,1:2] %>%
      tidyr::spread(characteristic_id,characteristic_value) %>%
      dplyr::mutate(comid=char_list$comid,
                    siteno=site)
    
    gage_chars <- rbind(gage_chars,temp_df)
  } else {
    print(paste("no data for", site))
  }
}

gage_chars <- mutate_at(gage_chars, vars(-comid,-siteno), funs(as.numeric))
write_feather(gage_chars,"data/basinchars/nldi/mobile_gage_chars.feather")

# extract basin characteristics for hucs
huc_chars <- NULL
for(i in 1:nrow(mobile)){
  print(i)
  huc <- mobile$HUC_12[i]
  char_list <- characteristics_nldi("huc12pp",huc)
  
  if(!is.null(char_list)) {
    temp_df <- char_list$characteristics[,1:2] %>%
      tidyr::spread(characteristic_id,characteristic_value) %>%
      dplyr::mutate(comid=char_list$comid,
                    HUC_12 = huc)
    
    huc_chars <- rbind(huc_chars,temp_df)
  } else {
    print(paste("no data for", huc))
  }
}

huc_chars <- mutate_at(huc_chars, vars(-comid,-HUC_12), funs(as.numeric))
write_feather(huc_chars,"data/basinchars/nldi/mobile_huc_chars.feather")

# DAYMET climate data ----
# create huc12 stencil for geoknife
huc_stencil <- as_Spatial(mobile$geometry,IDs=as.character(mobile$HUC_12)) %>%
  sp::spTransform(CRSobj = CRS("+proj=longlat +datum=WGS84")) %>%
  simplegeom()

# query available variables in prism dataset
daymet = webdata(url='https://thredds.daac.ornl.gov/thredds/dodsC/daymet-v3-agg/na.ncml')
query(daymet,'variables')

# load data
precip <- webdata(list(
  times = as.POSIXct(c('1980-01-01','2015-12-31')),
  url = 'http://cida-eros-netcdfdev.er.usgs.gov:8080/thredds/dodsC/thredds_workspace/daymet/daymet.ncml',
  variables = 'prcp'))

tmax <- webdata(list(
  times = as.POSIXct(c('1980-01-01','2015-12-31')),
  url = 'http://cida-eros-netcdfdev.er.usgs.gov:8080/thredds/dodsC/thredds_workspace/daymet/daymet.ncml',
  variables = 'tmax'))

tmin <- webdata(list(
  times = as.POSIXct(c('1980-01-01','2015-12-31')),
  url = 'http://cida-eros-netcdfdev.er.usgs.gov:8080/thredds/dodsC/thredds_workspace/daymet/daymet.ncml',
  variables = 'tmin'))

# precip <- webdata(list(
#   times = as.POSIXct(c('1980-01-01','2015-12-31')),
#   url = 'https://thredds.daac.ornl.gov/thredds/dodsC/daymet-v3-agg/na.ncml',
#   variables = 'prcp'))
# 
# tmax <- webdata(list(
#   times = as.POSIXct(c('1980-01-01','2015-12-31')),
#   url = 'https://thredds.daac.ornl.gov/thredds/dodsC/daymet-v3-agg/na.ncml',
#   variables = 'tmax'))
# 
# tmin <- webdata(list(
#   times = as.POSIXct(c('1980-01-01','2015-12-31')),
#   url = 'https://thredds.daac.ornl.gov/thredds/dodsC/daymet-v3-agg/na.ncml',
#   variables = 'tmin'))
# 
# vp <- webdata(list(
#   times = as.POSIXct(c('1980-01-01','2015-12-31')),
#   url = 'https://thredds.daac.ornl.gov/thredds/dodsC/daymet-v3-agg/na.ncml',
#   variables = 'vp'))
# 
# swr <- webdata(list(
#   times = as.POSIXct(c('1980-01-01','2015-12-31')),
#   url = 'https://thredds.daac.ornl.gov/thredds/dodsC/daymet-v3-agg/na.ncml',
#   variables = 'srad'))

# create webprocesses
precip_job <- geoknife(stencil=huc_stencil, fabric=precip, wait = F, email = 'scworland@usgs.gov')
tmax_job <- geoknife(stencil=huc_stencil, fabric=tmax, wait = F, email = 'scworland@usgs.gov')
tmin_job <- geoknife(stencil=huc_stencil, fabric=tmin, wait = F, email = 'scworland@usgs.gov')

# download data and save local copy
dmppt <- result(precip_job)
dmtmax <- result(tmax_job)
dmtmin <- result(tmin_job)
write_feather(dmppt,"data/dvs/daymet_ppt.feather")
write_feather(dmtmax,"data/dvs/daymet_tmax.feather")
write_feather(dmtmin,"data/dvs/daymet_tmin.feather")






