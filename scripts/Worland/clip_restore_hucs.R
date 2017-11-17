
library(sf)

# Download WBD data here: http://www.horizon-systems.com/NHDPlus/V2NationalData.php

# load full HU12s dataset
huc12 <- st_read("data/nhd_wbd/WBDSnapshot_National.shp", stringsAsFactors = F)

# load clipping extent
bounds <- st_read("data/restore_hucs/extent/extent4.shp", stringsAsFactors = F)

# set CRS equal
st_crs(bounds) <- st_crs(huc12)

# find intersection
intersect_restore <- st_intersection(huc12, bounds) 

# write to shape file
st_write(intersect_restore,"data/restore_hucs/restoreHUC12s.shp")


