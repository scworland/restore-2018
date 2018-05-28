
library(pacman)
pacman::p_load(tidyverse,stringr,lubridate,feather,geoknife,sf,sp,dataRetrieval,ggsn)
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

# site map
pplo <- read_csv("data/huc12/spCOV.csv") %>%
  select(lat,lon,est_pplo) %>%
  group_by(lat,lon) %>%
  summarize(pplo = mean(est_pplo,na.rm=T)) 

boundary <- read_sf("data/shapefiles/simple_bounds/simple_bounds.shp")
gages <- read_sf("data/shapefiles/gages/gage_list_comid.shp")
hucpp <- read_feather("data/huc12/all_huc12_covariates.feather") %>%
  select(lat,lon,runoff_mean) %>%
  group_by(lat,lon) %>%
  summarize(Q = mean(runoff_mean)) 

regions = c("texas","alabama","florida","mississippi",
            "louisiana","georgia","arkansas","tennessee",
            "oklahoma","missouri","kentucky")

states <- subset(map_data("state"),region %in% regions)

ggplot() +
  geom_polygon(data=states, aes(x=long, y=lat, group=group),fill="grey85",
               linetype="dotted",color="grey20",size=0.3) +
  geom_sf(data=boundary,fill="white",alpha=0.7,color="grey15") +
  #geom_point(data=pplo,aes(x=lon,y=lat,color=pplo),size=0.1,show.legend=T) +
  #scale_color_viridis_c(option="plasma",name="p(Q=0)") +
  #geom_sf(data=gages,shape=4) +
  coord_sf(crs = st_crs(gages), datum = NA) +
  north(states, symbol = 3, scale = 0.10, location="topleft") +
  scalebar(states, dist = 250, dd2km = TRUE, model = 'WGS84', 
           st.size = 3, location="bottomleft") +
  theme_void() 



