library(rgeos)
library(leaflet)
library(rgdal)

source("scripts/worland/nldi_funs.R")

query_nldi("ss_gages","08116650",tier = "test")
query_nldi("nwissite","USGS-03568400")

nldiURLs <- list(site_data = "https://cida-test.er.usgs.gov/nldi/ss_gages/08116650",
                 basin_boundary = "https://cida-test.er.usgs.gov/nldi/ss_gages/08116650/basin",
                 UT = "https://cida-test.er.usgs.gov/nldi/ss_gages/08116650/navigate/UT",
                 UM = "https://cida-test.er.usgs.gov/nldi/ss_gages/08116650/navigate/UM",
                 DM = "https://cida-test.er.usgs.gov/nldi/ss_gages/08116650/navigate/DM")

nldi_data <- list()

for(n in names(nldiURLs)) {
  nldi_data[n] <- rgdal::readOGR(dsn = nldiURLs[n][[1]], layer = "OGRGeoJSON", verbose = FALSE)
  print(paste(n, "is of class", class(nldi_data[n][[1]])[1], "and has", length(nldi_data[n][[1]]), "features"))
}


map <- get_base_map() # the function described above.

map <- leaflet::addPolygons(map, 
                            data=nldi_data$basin_boundary, 
                            color = "orange", 
                            fill = FALSE, 
                            weight = 2,
                            opacity = 1)

map <- leaflet::addPolylines(map, 
                             data = nldi_data$UT,
                             color = "blue",
                             weight = 1,
                             opacity = 1)

map <- leaflet::addPolylines(map,
                             data = nldi_data$UM, 
                             color = "blue", 
                             weight = 3, 
                             opacity = 0.5)

map <- leaflet::addCircleMarkers(map = map,
                                 data = nldi_data$UTwqp,
                                 radius = 1,
                                 color = "black",
                                 opacity = .5,
                                 fill = FALSE, 
                                 popup = UTwqp_html)

map <- leaflet::addCircleMarkers(map, 
                                 data = nldi_data$site_data, 
                                 radius = 5, 
                                 color = "red")

map <- leaflet::addPolylines(map,
                             data = nldi_data$DM, 
                             color = "blue", 
                             weight = 3, 
                             opacity = 0.5)
