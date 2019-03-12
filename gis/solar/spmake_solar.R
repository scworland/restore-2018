library(sp)
library(rgdal)

ALBEA <- paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 ",
                "+x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
ALBEA <- sp::CRS(ALBEA)

unzip("us9809_dni_updated.zip", overwrite=TRUE)
spDNI_1998to2009 <- readOGR("us9809_dni_updated", "us9809_dni_updated")
unlink("./us9809_dni_updated/", recursive=TRUE)
if(dir.exists("./__MACOSX")) unlink("./__MACOSX", recursive=TRUE)
save(spDNI_1998to2009, file="spDNI_1998to2009.RData")

unzip("lower_48_dni_10km.zip", overwrite=TRUE)
spDNI_1998to2005 <- readOGR("lower_48_dni_10km", "us9805_dni")
unlink("./lower_48_dni_10km/", recursive=TRUE)
if(dir.exists("./__MACOSX")) unlink("./__MACOSX", recursive=TRUE)
save(spDNI_1998to2005, file="spDNI_1998to2005.RData")


