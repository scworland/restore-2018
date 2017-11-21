# https://cida.usgs.gov/nldi/huc12pp/170900120202/basin
#' @example basin <- get_basin_nldi("huc12pp", "170900120202")
get_basin_nldi <- function(f_source, f_id, tier = "prod") {
  nldi_base_url <- "https://cida.usgs.gov/nldi"
  
  url <- paste(nldi_base_url, f_source, f_id, "basin", 
               sep = "/")
  
  c <- ""
  
  try(c <- rawToChar(httr::GET(url)$content), silent = T)
  
  if(nchar(c)==0) {
    NULL
  } else {
    sf::st_read(c, quiet = T, stringsAsFactors = F)
  }
}

hucs <- c("070200121110", "070700051701", "031601130201", "160201020603")

basins <- setNames(as.list(1:length(hucs)), hucs)

for(huc in hucs) {
  basins[huc] <- get_basin_nldi("huc12pp", huc)
}

# I know there is a better way...
# Just need to call c() with all elements of the list in there.
basins_sf <- basins[[1]]
for(i in 2:length(hucs)) {
  basins_sf <- c(basins_sf, basins[[i]])  
}

# Can create an sf object with the huc12s as attributes like:
basins_sf <- sf::st_sf(geometry = basins_sf, huc12 = hucs, stringsAsFactors = F, row.names = hucs)

# Then convert to a geoknife ready sp object:
basins_sp <- sp::spTransform(sp::SpatialPolygonsDataFrame(Sr = sf::as_Spatial(basins_sf$geometry), 
                                                          data = data.frame(huc12=basins_sf$huc12, stringsAsFactors = F), match.ID = F), CRSobj = CRS("+proj=longlat +datum=WGS84"))

# Now use geoknife:
stencil <- geoknife::simplegeom(basins_sp)


