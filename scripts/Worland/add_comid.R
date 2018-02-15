library(pacman)
pacman::p_load(tidyverse,stringr,lubridate,feather,geoknife,sf,sp,dataRetrieval)

# uncomment to run from command line
# setwd("/Users/scworlan/Documents/Restore")

source('scripts/worland/nldi_funs.R')
source('scripts/worland/utils.R')

# load HUC12s for restore footprint
rest_huc12 <- st_read('data/shapefiles/restore_hucs/restoreHUC12s.shp', stringsAsFactors = F)

# load site list and get NHD comids for each gage
sites <- read_csv("data/decade1950plus_site_list.csv")

UM = "https://cida-test.er.usgs.gov/nldi/ss_gages/08116650/navigate/UM"
nldi_data <- rgdal::readOGR(dsn = UM, layer = "OGRGeoJSON", verbose = FALSE)

# use NLDI to find comids
sites$comid <- NA
for(i in 1:nrow(sites)){
  d <- query_nldi("ss_gages",sites$site_no[i],tier = "test")
  if(!is.null(d)) {
    sites$comid[i] <- d$features$properties$comid[1]
  } else {
    print(paste0("no data for ", sites$site_no[i]))
  }
  print(paste0("site #: ", i))
}

sites <- sites %>% 
  na.omit()

write_csv(sites,"data/decade1950plus_site_list_comid.csv")
write_feather(sites,"data/decade1950plus_site_list.feather")

# ----test COMIDs-----
mwcomid <- read.dbf("data/gage/Gages_Order_PathIDs.dbf", as.is = FALSE)
cda <- read_delim("data/gage/cda_bust.txt",delim=" ") %>%
  mutate(bust="cda_bust")

sites <- read_feather("data/decade1950plus_site_list.feather") %>%
  left_join(mwcomid,by="site_no") %>%
  select(site_no,nldi_comid=comid,gageloc_comid=GL_COMID) %>%
  mutate_all(funs(as.character)) %>%
  mutate(diff=ifelse(nldi_comid!=gageloc_comid,"1","0"),
         diff=replace(diff,is.na(diff),"1")) %>%
  filter(diff==1) %>%
  left_join(cda,by="site_no")

write_feather(sites,"data/gage/comid_diff.feather")

# find coordinates of gages 
# sites <- read_feather("data/decade1950plus_site_list.feather")

site_info <- readNWISdata(sites=sites$site_no, service="site") %>%
  select(site_no,lon=dec_long_va,lat=dec_lat_va) %>%
  left_join(sites, by="site_no") %>%
  st_as_sf(coords=c("lon","lat"), crs=st_crs(rest_huc12), remove=F) %>%
  st_intersection(rest_huc12) %>%
  select(site_no,comid,lon,lat,HUC_12) %>%
  distinct(site_no,.keep_all=T) %>%
  na.omit() 

# get NHD comids for each huc 12
huc12_list <- data.frame(huc12=rest_huc12$HUC_12) %>%
  mutate(huc12 = as.character(huc12))

# use NLDI to find comids
huc12_list$comid <- NA
for(i in 1:nrow(huc12_list)){
  d <- query_nldi("huc12pp",huc12_list$huc12[i],tier = "test")
  if(!is.null(d)) {
    huc12_list$comid[i] <- d$features$properties$comid[1]
  } else {
    print(paste0("no data for ", huc12_list$huc12[i]))
  }
  print(paste0("huc #: ", i))
}

huc12_list <- huc12_list %>%
  mutate(comid = ifelse(comid == "",NA,comid)) %>%
  na.omit()

write_csv(huc12_list,"data/huc12_list_comid.csv")
write_feather(huc12_list,"data/huc12_list_comid.feather")

# generate shapefiles
st_write(site_info, 'data/shapefiles/gages/restore_gages.shp')
st_write(rest_huc12, 'data/shapefiles/restore_hucs/restoreHUC12s.shp')


# # create immutable basin characteristic matrix
# items <- c("5783e986e4b0ac0b97f5428c",
#            "57976a0ce4b021cadec97890",
#            "57c5f12fe4b0f2f0cebdab0d")
# 
# immutable_chars <- list()
# for (i in 1:length(items)){
#   immutable_chars[[i]] <- sw_sb_extract(items[i])
# }
# 
# 
# hold <- Reduce(merge,immutable_chars)

