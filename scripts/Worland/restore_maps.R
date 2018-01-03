
library(pacman)
pacman::p_load(tidyverse,stringr,lubridate,feather,geoknife,sf,sp,dataRetrieval)
source('scripts/worland/nldi_funs.R')

# NHD+ basin characteristics ----
# load HUC12s for restore footprint
rest_huc12 <- st_read('data/shapefiles/restore_hucs/restoreHUC12s.shp', stringsAsFactors = F)

# subset mobile basin
mobile <- filter(rest_huc12, str_detect(HUC_8, '0316')) %>%
  distinct(HUC_12, .keep_all=T)

basin <- get_basin_nldi("huc12pp", mobile$HUC_12[8])

ggplot() + 
  geom_sf(data=basin,aes(),fill=NA,linetype="dashed") + 
  geom_sf(data=mobile[8,], aes(),fill="blue",alpha=0.5,color=NA) +
  theme_bw() +
  ggtitle("HUC 031601080208")

