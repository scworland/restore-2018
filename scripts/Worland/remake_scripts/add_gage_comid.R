
add_gage_comid <- function(site_list,gage_locs,updated_locs,wbd_huc12){
  
  # subset columns and make df
  gage_locs <- gage_locs %>% 
    select(site_no=SOURCE_FEA,comid=FLComID) %>%
    st_set_geometry(NULL)
  
  # make df
  updated_locs <- st_set_geometry(updated_locs,NULL)
  
  # join to site list
  sites_comids <- site_list %>%
    data.frame() %>%
    left_join(gage_locs,by="site_no") %>%
    left_join(updated_locs,by=c("site_no"="Gage_no")) %>%
    mutate(COMID = ifelse(COMID=="-9999",NA,COMID),
           comid = ifelse(!is.na(COMID),COMID,comid)) %>%
    select(site_no,comid) %>%
    mutate(comid=as.character(comid)) 
  
  # use NWIS to get more site information
  gage_list_comid <- readNWISdata(sites=sites_comids$site_no, service="site") %>%
    select(site_no,lon=dec_long_va,lat=dec_lat_va) %>%
    left_join(sites_comids, by="site_no") %>%
    st_as_sf(coords=c("lon","lat"), crs=st_crs(wbd_huc12), remove=F) %>%
    st_intersection(wbd_huc12) %>% 
    select(site_no,comid,huc12=HUC_12,lon,lat) %>%
    distinct(site_no,.keep_all=T) %>%
    st_set_geometry(NULL)
  
  return(gage_list_comid)
}