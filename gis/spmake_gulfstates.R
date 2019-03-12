library(sp)
library(rgdal)

unzip("GulfStates.zip", overwrite=TRUE)
GulfStates          <- readOGR("GulfStates/", "GulfStates")
GulfStates_modified <- readOGR("GulfStates/", "GulfStates_modified")
unlink("./GulfStates/", recursive=TRUE)
if(dir.exists("./__MACOSX")) unlink("./__MACOSX", recursive=TRUE)

ALBEA <- paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 ",
     "+x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
ALBEA <- sp::CRS(ALBEA)

GulfStates          <- spTransform(GulfStates,          ALBEA)
GulfStates_modified <- spTransform(GulfStates_modified, ALBEA)

save(GulfStates, GulfStates_modified, file="GulfStates.RData")
