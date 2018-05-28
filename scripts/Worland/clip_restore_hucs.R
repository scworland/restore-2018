
library(sf)

# Download WBD data here: http://www.horizon-systems.com/NHDPlus/V2NationalData.php

# load full HU12s dataset
huc12 <- st_read("data/shapefiles/nhd_wbd/WBDSnapshot_National.shp", stringsAsFactors = F)

# load clipping extent
bounds <- st_read("data/shapefiles/restore_bnd/restore_bnd.shp", stringsAsFactors = F) 

# set CRS equal
st_crs(bounds) <- st_crs(huc12)

# find intersection
restore_hucs <- st_intersection(huc12,bounds)

# write to shape file
st_write(restore_hucs,"data/shapefiles/restore_hucs/restoreHUC12s.shp")


