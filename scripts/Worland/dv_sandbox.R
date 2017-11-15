
library(tidyverse)
library(dataRetrieval)
library(lubridate)
library(feather)

# laod R environment
load("data/dvs/DV.RData")

# extract site numbers
sites <- ls(DV)

# load list of dvs
dv_list <- as.list(DV)

# load site descriptions (lat, lon, area, period of record)
site_desc <- read_feather("data/basinchars/BASIN_CHAR_TOT.feather") %>%
  mutate(siteno=sprintf("%08d",siteno),
         min_year = year(min_date),
         max_year = year(max_date)) %>%
  select(siteno,da_km2=DASqKm,lat=dec_lat_va,
         lon=dec_long_v,min_year,max_year)

# map sites
state <- map_data('state') 

ggplot() + geom_polygon(data=state,aes(long,lat,group=group),color = "black", fill= "white",size=0.5) +
  coord_fixed(1.3) + theme_bw() + xlab("") + ylab("") + 
  geom_point(data=site_desc,aes(lon,lat),fill="dodgerblue",shape=21) +
  theme_void() + 
  coord_fixed(1.3,xlim=c(-106,-75),ylim=c(23,40))
