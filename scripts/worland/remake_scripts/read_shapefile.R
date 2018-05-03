 
sw_shp_read <- function(file) {
  out <- st_read(file, stringsAsFactors = F)
  return(out)
}