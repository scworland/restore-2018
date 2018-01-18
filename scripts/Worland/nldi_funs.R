query_nldi <- function(f_source, f_id, tier = "prod") {
  nldi_base_url <- get_nldi_url(tier)
  
  url <- paste(nldi_base_url, f_source, f_id, 
               sep = "/")
  
  c <- rawToChar(httr::GET(url)$content)
  
  if(nchar(c)==0) {
    NULL
  } else {
    try(jsonlite::fromJSON(c), silent = F)
  }
}

navigate_nldi <- function(f_source, f_id, mode = "UM", 
                          d_source = NULL, distance = NULL, tier = "prod") {
  nldi_base_url <- get_nldi_url(tier)
  
  url <- paste(nldi_base_url, f_source, f_id, "navigate", mode, d_source, 
               sep = "/")
  
  if(!is.null(distance)) {
    url <- paste0(url, "?distance=", distance)
  }
  
  c <- rawToChar(httr::GET(url)$content)
  
  if(nchar(c)==0) {
    NULL
  } else {
    try(jsonlite::fromJSON(c), silent = F)
  }
  
}

sapply_char_fun <- function(x) as.numeric(characteristics_nldi(gf_source_id, x, char_type, char_id)$characteristics$characteristic_value)

# https://cida.usgs.gov/nldi/huc12pp/070900020503/TOT?characteristicId=TOT_BASIN_AREA
characteristics_nldi <- function(f_source, f_id, char_type = "TOT", char_id = NULL, tier = "prod") {
  nldi_base_url <- get_nldi_url(tier)
  
  url <- paste(nldi_base_url, f_source, f_id, char_type, 
               sep = "/")
  
  if(!is.null(char_id)) {
    url <- paste0(url, "?characteristicId=", char_id)
  }
  
  c <- ""
  r <- httr::GET(url)
  
  if(r$status_code != 400) {
    try(c <- rawToChar(r$content), silent = T)
  }
  
  if(nchar(c)==0) {
    NULL
  } else {
    try(jsonlite::fromJSON(c), silent = F)
  }
}

# https://cida.usgs.gov/nldi/huc12pp/170900120202/basin
#' @example basin <- get_basin_nldi("huc12pp", "170900120202")
get_basin_nldi <- function(f_source, f_id, tier = "prod") {
  
  nldi_base_url <- get_nldi_url(tier)
  
  url <- paste(nldi_base_url, f_source, f_id, "basin", 
               sep = "/")
  
  c <- ""
  
  try(c <- rawToChar(httr::GET(url)$content), silent = T)
  
  if(nchar(c)==0) {
    NULL
  } else {
    try(sf::st_read(c, quiet = T, stringsAsFactors = F), silent = F)
  }
}

get_nldi_url <- function(tier = "prod") {
  if(tier=="prod") {
      "https://cida.usgs.gov/nldi"
  } else if(tier=="test") {
    "https://cida-test.er.usgs.gov/nldi"
  } else if(tier=="local") {
    "http://localhost:8080/nldi"
  }
}


GetURL <- function(service, host = "basemap.nationalmap.gov") {
  sprintf("https://%s/arcgis/services/%s/MapServer/WmsServer", host, service)
}

get_base_map <- function(options = leaflet::leafletOptions()) {
  map <- leaflet::leaflet(options = options)
  grp <- c("USGS Topo", "USGS Imagery Only", "USGS Imagery Topo",
           "USGS Shaded Relief", "Hydrography")
  att <- paste0("<a href='https://www.usgs.gov/'>",
                "U.S. Geological Survey</a> | ",
                "<a href='https://www.usgs.gov/laws/policies_notices.html'>",
                "Policies</a>")
  map <- leaflet::addWMSTiles(map, GetURL("USGSTopo"),
                              group = grp[1], attribution = att, layers = "0")
  map <- leaflet::addWMSTiles(map, GetURL("USGSImageryOnly"),
                              group = grp[2], attribution = att, layers = "0")
  map <- leaflet::addWMSTiles(map, GetURL("USGSImageryTopo"),
                              group = grp[3], attribution = att, layers = "0")
  map <- leaflet::addWMSTiles(map, GetURL("USGSShadedReliefOnly"),
                              group = grp[4], attribution = att, layers = "0")
  opt <- leaflet::WMSTileOptions(format = "image/png", transparent = TRUE)
  map <- leaflet::addWMSTiles(map, GetURL("USGSHydroCached"),
                              group = grp[5], options = opt, layers = "0")
  map <- leaflet::hideGroup(map, grp[5])
  opt <- leaflet::layersControlOptions(collapsed = FALSE)
  map <- leaflet::addLayersControl(map, baseGroups = grp[1:4],
                                   overlayGroups = grp[5], options = opt)
}
