
library(rgdal)
library(sp)

# Download WBD data here: http://www.horizon-systems.com/NHDPlus/V2NationalData.php
# Read in full HU12s dataset
huc12 <- readOGR(path.expand("data/nhd_wbd"), "WBDSnapshot_National")

# Read in table of HU12s in restore footprint
restore12 <- base::subset(huc12, HUC_12 %in% c("GA","AL,GA","FL,GA","AL,FL,GA"))

