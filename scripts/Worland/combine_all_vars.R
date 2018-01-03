
library(pacman)
pacman::p_load(tidyverse,stringr,lubridate,feather,sf)
source("scripts/worland/utils.R")

# Create gages dataset ----

## load immutable variables
immutable_gage_df <- read_feather("data/basinchars/nhd_sb/immutable_gage.feather")

## load mutable variables
climate_gage_df <- read_feather("data/basinchars/nhd_sb/climate_gage_df.feather")
house_gage_df <- read_feather("data/basinchars/nhd_sb/housedens_gage_df.feather")
dam_gage_df <- read_feather("data/basinchars/nhd_sb/dam_gage_df.feather")

lulc_gage_df <- read_feather("data/basinchars/nhd_sb/lulc_gage_df.feather") %>%
  filter(decade != 1940)

hold <- list(climate_gage_df,house_gage_df,dam_gage_df,lulc_gage_df) %>%
  reduce(left_join, by = c("comid","site_no","decade")) %>%
  arrange(comid,decade) %>%
  left_join(immutable_gage_df, by = c("comid","site_no"))