
library(pacman)
pacman::p_load(tidyverse,stringr,lubridate,feather,sf)

# Load covariate data
hucs <- st_read('data/shapefiles/mobile_shpfiles/mobileHUC12s.shp', stringsAsFactors = F)
sites <- st_read('data/shapefiles/mobile_shpfiles/mobile_gages.shp', stringsAsFactors = F)
huc_chars <- read_feather("data/basinchars/nldi/mobile_huc_chars.feather")
gage_chars <- read_feather("data/basinchars/nldi/mobile_gage_chars.feather")
dmppt <- read_feather("data/dvs/daymet_ppt.feather")
dmtmax <- read_feather("data/dvs/daymet_tmax.feather")
dmtmin <- read_feather("data/dvs/daymet_tmin.feather")

# Load daily streamflow



hucs2 <- hucs %>%
  filter(AreaHUC12 < max(AreaHUC12)) %>%
  left_join(huc_chars, by="HUC_12")

plot(hucs2[c("AreaHUC12","TOT_BASIN_AREA","CAT_AREA_SQKM")])

plot(hucs2[c("TOT_ET","TOT_ELEV_MEAN","TOT_RECHG")])
