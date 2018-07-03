

add_huc12_comid_pp <- function(restore_hucs,huc12_pp) {
  
  # transform huc 12 pp
  huc12_pp <- st_transform(huc12_pp,"+proj=longlat +datum=WGS84")
  
  # add coordinates as columns
  pp_coords <- st_coordinates(huc12_pp)
  huc12_pp$lon <- pp_coords[,1]
  huc12_pp$lat <- pp_coords[,2]
  
  # turn off pour point geometry
  st_geometry(huc12_pp) <- NULL
  
  # add pour points to huc12 and write file
  huc12_comid_pp <- restore_hucs %>%
    select(HUC_12) %>%
    left_join(huc12_pp,by="HUC_12") %>%
    mutate(comid=as.character(COMID)) %>%
    select(huc12=HUC_12,comid,lon,lat)
  
  # turn off huc12 polygon geometry 
  st_geometry(huc12_comid_pp) <- NULL
  
  return(huc12_comid_pp)
}