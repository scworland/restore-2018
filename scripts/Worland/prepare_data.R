
library(pacman)
pacman::p_load(tidyverse,stringr,lubridate,feather,sf)
source("scripts/worland/utils.R")

# Load covariate data
hucs <- st_read('data/shapefiles/mobile_shpfiles/mobileHUC12s.shp', stringsAsFactors = F)
site_info <- st_read('data/shapefiles/mobile_shpfiles/mobile_gages.shp', stringsAsFactors = F)
huc_chars <- read_feather("data/basinchars/nldi/mobile_huc_chars.feather")
gage_chars <- read_feather("data/basinchars/nldi/mobile_gage_chars.feather")
dmppt <- read_feather("data/dvs/daymet_ppt.feather")
dmtmax <- read_feather("data/dvs/daymet_tmax.feather")
dmtmin <- read_feather("data/dvs/daymet_tmin.feather")

# change shape of climate data 
ppt <- dmppt %>%
  rename(date=DateTime) %>%
  mutate(date = as.Date(date, "%d/%m/%y",tz="UTC")) %>%
  gather(HUC_12, ppt, -date) 

tmax <- dmtmax %>%
  rename(date=DateTime) %>%
  mutate(date = as.Date(date, "%d/%m/%y",tz="UTC")) %>%
  gather(HUC_12, tmax, -date) 

tmin <- dmtmin %>%
  rename(date=DateTime) %>%
  mutate(date = as.Date(date, "%d/%m/%y", tz="UTC")) %>%
  gather(HUC_12, tmin, -date) 

# build model data using gages ----

# Load daily streamflow
load("data/dvs/DV.RData")
sites <- ls(DV)
dv_list <- as.list(DV)
dv_list <- remove_aux(dv_list) 

# subset for mobile gages and daymet years
mobile_dvs <- dv_list[sites[which(sites %in% site_info$siteno)]]

dvs <- do.call("rbind", mobile_dvs) %>%
  filter(year >= 1980 & year <= 2015) %>%
  select(siteno=site_no,date=Date,Q=Flow,
         cd=Flow_cd,year,month,decade,wyear)

gage_time <- dvs %>%
  left_join(site_info,by="siteno") %>%
  left_join(ppt,by=c("HUC_12","date")) %>%
  left_join(tmin,by=c("HUC_12","date")) %>%
  left_join(tmax,by=c("HUC_12","date")) %>%
  select(siteno,lat,lon,HUC_12,date,
         year,month,Q,ppt,tmin,tmax)

gage_static <- gage_chars %>%
  rename_all(tolower) %>%
  select(siteno,comid,everything())

write_feather(gage_time,"data/gage/gage_time.feather")
write_feather(gage_static,"data/gage/gage_static.feather")

# build prediction data using hucs ----
huc12_time <- ppt %>%
  left_join(tmin,by=c("HUC_12","date")) %>%
  left_join(tmax,by=c("HUC_12","date")) %>%
  mutate(year = year(date),
         month = month(date)) %>%
  select(HUC_12,date,year,month,ppt,tmin,tmax)

huc12_static <- huc_chars %>%
  rename_at(vars(-HUC_12),tolower) %>%
  select(HUC_12,comid,everything())

write_feather(huc12_time,"data/huc12/huc12_time.feather")
write_feather(huc12_static,"data/huc12/huc12_static.feather")


