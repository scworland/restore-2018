
clip_restore_hucs <- function(restore_boundary,wbd_huc12) {
  # set CRS 
  st_crs(restore_boundary) <- st_crs(wbd_huc12)
  
  # find intersection
  restore_hucs <- st_intersection(wbd_huc12,restore_boundary)
  
  return(restore_hucs)
}




