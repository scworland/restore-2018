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

