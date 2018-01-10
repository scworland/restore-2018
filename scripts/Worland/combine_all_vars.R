
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

lulc_gage_df <- read_feather("data/basinchars/nhd_sb/lulc_gage_df.feather") 

all_gage_covariates <- list(climate_gage_df,house_gage_df,dam_gage_df,lulc_gage_df) %>%
  reduce(left_join, by = c("comid","site_no","decade")) %>%
  arrange(comid,decade) %>%
  left_join(immutable_gage_df, by = c("comid","site_no")) %>%
  rename_all(tolower) %>%
  select(comid:woody_wetland,tot_bfi:tot_twi,tot_basin_area:tot_rdx,everything()) %>%
  filter(decade %in% c(1950,1960,1970,1980,1990,2000))

## load FDC data from Will A. 
fdc50_00 <- read_feather("data/gage/fdc_lmr_pplo.feather")  %>%
  rename(site_no=site)

## combine everything
all_gage_data <- fdc50_00 %>%
  left_join(all_gage_covariates, by = c("site_no","decade")) %>%
  select(comid,site_no,decade,everything()) %>%
  arrange(site_no,decade) %>%
  mutate(decade = as.character(decade)) %>%
  filter(decade != 2010) %>%
  na.omit()

write_feather(all_gage_data,"data/gage/all_gage_data.feather")

# Create huc12 dataset ----

## load immutable variables
immutable_huc12_df <- read_feather("data/basinchars/nhd_sb/immutable_huc12.feather")

## load mutable variables
climate_huc12_df <- read_feather("data/basinchars/nhd_sb/climate_huc12_df.feather")
house_huc12_df <- read_feather("data/basinchars/nhd_sb/housedens_huc12_df.feather")
dam_huc12_df <- read_feather("data/basinchars/nhd_sb/dam_huc12_df.feather")

lulc_huc12_df <- read_feather("data/basinchars/nhd_sb/lulc_huc12_df.feather") 

all_huc12_covariates <- list(climate_huc12_df,house_huc12_df,dam_huc12_df,lulc_huc12_df) %>%
  reduce(left_join, by = c("comid","huc12","decade")) %>%
  arrange(comid,decade) %>%
  left_join(immutable_huc12_df, by = c("comid","huc12")) %>%
  filter(decade %in% c(1950,1960,1970,1980,1990,2000))


write_feather(all_huc12_covariates,"data/huc12/all_huc12_covariates.feather")

