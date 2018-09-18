
clip_restore_hucs <- function(restore_pp,wbd_huc12) {
  
  # find huc12 polygons in boundary
  restore_hucs <- wbd_huc12 %>%
    filter(HUC_12 %in% restore_pp$HUC_12)
  
  return(restore_hucs)
}




