
library(tidyverse)
library(dataRetrieval)
library(lubridate)
library(feather)
source("scripts/worland/utils.R")

# laod R environment
load("data/dvs/DV.RData")

# extract site numbers
sites <- ls(DV)

# load list of dvs
dv_list <- as.list(DV)

# LULC files
pth = "data/basinchars"
lulc_files <- list.files(path=pth, pattern = glob2rx("*LULC*.feather"),full.names = T)
temp_files <- list.files(path=pth, pattern = glob2rx("*TAV*.feather"),full.names = T)[2:13]
ppt_files <- list.files(path=pth, pattern = glob2rx("*PPT*.feather"),full.names = T)

temp <- read_months(temp_files,"temp") %>%
  reduce(left_join, by = "siteno")

ppt <- read_months(ppt_files,"ppt") %>%
  reduce(left_join, by = "siteno")

lulc <- read_years(lulc_files) %>%
  reduce(left_join, by = "siteno") %>%
  select(matches('Cult|Dev|Dec|Ever|Hay|Mixed')) %>%
  select(contains('2004'))
  
# load basin chars
bfi <- read_feather("data/basinchars/BFI.feather")
et <- read_feather("data/basinchars/ET.feather")
runoff <- read_feather("data/basinchars/RUN7100.feather")
tav <- read_feather("data/basinchars/TAV7100_ANN.feather")

df_list <- list(bfi,et,runoff,tav) 

X <- df_list %>%
  reduce(left_join, by = "siteno") %>%
  select(siteno,contains("CAT")) %>%
  left_join(temp) %>% 
  left_join(ppt) %>%
  left_join(lulc)

#hold <- read.dbf("data/basinchars/geodatabases/lulc_2003.dbf")

# load site descriptions (lat, lon, area, period of record)
site_desc <- read_feather("data/basinchars/BASIN_CHAR_TOT.feather") %>%
  mutate(siteno=sprintf("%08d",siteno),
         min_year = year(min_date),
         max_year = year(max_date)) %>%
  select(siteno,da_km2=DASqKm,lat=dec_lat_va,
         lon=dec_long_v,min_year,max_year)


# map sites
states <- map_data('state', region=c("texas","georgia","florida",
                                     "louisiana","alabama", "tennessee",
                                     "mississippi","arkansas","oklahoma",
                                     "kentucky","missouri")) 

ggplot() + geom_polygon(data=states,aes(long,lat,group=group),color = "black", fill= "white",size=0.5) + 
  geom_point(data=site_desc,aes(lon,lat),fill="dodgerblue",shape=21,alpha=0.5) +
  coord_fixed(1.3) + theme_void()
