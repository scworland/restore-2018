library(sp)
library(rgdal)

unzip("restore_huc12.zip", overwrite=TRUE)
huc12 <- readOGR("restore_huc12/", "restore_huc12")
huc12_projection <- proj4string(huc12)
unlink("./restore_huc12/", recursive=TRUE)
if(dir.exists("./__MACOSX")) unlink("./__MACOSX", recursive=TRUE)

ALBEA <- paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 ",
     "+x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
ALBEA <- sp::CRS(ALBEA)

huc12 <- spTransform(huc12, ALBEA)

save(huc12, file="restore_huc12.RData")
