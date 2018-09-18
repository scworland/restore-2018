

add_huc12_comid_pp <- function(restore_pp) {
  
  # subset coords, comid, and huc12
  huc12_comid_pp <- restore_pp %>%
    mutate(lon = st_coordinates(.)[,1],
           lat = st_coordinates(.)[,2],
           comid = as.character(COMID)) %>%
    select(huc12=HUC_12,comid,lon,lat) %>%
    st_set_geometry(NULL)
  
  return(huc12_comid_pp)
}