
sw_clip_restore_pp <- function(restore_boundary,huc12_pp) {
  
  # set CRS 
  huc12_pp <- st_transform(huc12_pp,crs=st_crs(restore_boundary))
  
  # find pour points within boundary
  pp_id <- st_within(huc12_pp,restore_boundary) %>%
    as.data.frame()
  
  restore_pp <- huc12_pp[pp_id$row.id,]
  
  return(restore_pp)
}