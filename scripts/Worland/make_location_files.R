
library(sf)
library(tidyverse)
library(feather)
library(dataRetrieval)

# load HUC12s for restore footprint
rest_huc12 <- st_read("data/shapefiles/restore_hucs/restoreHUC12s.shp", stringsAsFactors = F)

# load huc12 pour points
huc12_pp <- st_read("data/shapefiles/huc_pp_comid/HUC12_Outlet_COMIDs_CONUS.shp", stringsAsFactors = F) %>%
  st_transform("+proj=longlat +datum=WGS84")

# add coordinates as columns
pp_coords <- st_coordinates(huc12_pp)
huc12_pp$lon <- pp_coords[,1]
huc12_pp$lat <- pp_coords[,2]

# turn off pour point geometry
st_geometry(huc12_pp) <- NULL

# add pour points to huc12 and write file
huc12_list_comid <- rest_huc12 %>%
  left_join(huc12_pp,by="HUC_12") %>%
  mutate(COMID=as.character(COMID)) %>%
  select(huc12=HUC_12,comid=COMID,lon,lat)

st_write(huc12_list_comid,"data/shapefiles/restore_huc_comid/huc12_list_comid.shp")

# turn off huc12 polygon geometry and write feather
st_geometry(huc12_list_comid) <- NULL

write_feather(huc12_list_comid,"data/huc12_list_comid.feather")

# load site list and get NHD comids for each gage
gage_locs <- st_read("data/shapefiles/GageLoc/GageLoc.shp") %>% 
  select(site_no=SOURCE_FEA,comid=FLComID)
st_geometry(gage_locs) <- NULL

updated_locs <- st_read("data/shapefiles/CATCHMENT_gageloc_v1/CATCHMENT_gageloc_v1.shp")
st_geometry(updated_locs) <- NULL

# join to site list
sites <- read_csv("data/decade1950plus_site_list.csv") %>%
  data.frame() %>%
  left_join(gage_locs,by="site_no") %>%
  left_join(updated_locs,by=c("site_no"="Gage_no")) %>%
  mutate(comid = ifelse(!is.na(COMID),COMID,comid)) %>%
  select(site_no,comid) %>%
  filter(comid != -9999) %>%
  mutate(comid=as.character(comid)) %>%
  na.omit()

# use NWIS to get more site information
gage_list_comid <- readNWISdata(sites=sites$site_no, service="site") %>%
  select(site_no,lon=dec_long_va,lat=dec_lat_va) %>%
  left_join(sites, by="site_no") %>%
  st_as_sf(coords=c("lon","lat"), crs=st_crs(rest_huc12), remove=F) %>%
  st_intersection(rest_huc12) %>%
  select(site_no,huc12=HUC_12,comid,lon,lat) %>%
  distinct(site_no,.keep_all=T) %>%
  na.omit() 

# write shapefile
st_write(gage_list_comid,"data/shapefiles/gages/gage_list_comid.shp")

st_geometry(gage_list_comid) <- NULL

write_feather(gage_list_comid,"data/gage_list_comid.feather")







