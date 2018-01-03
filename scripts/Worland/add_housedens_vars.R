
library(pacman)
pacman::p_load(tidyverse,stringr,lubridate,feather,geoknife,sf,sp,dataRetrieval)

# uncomment to run from command line
setwd("/Users/scworlan/Documents/Restore")

# source scripts
source('scripts/worland/nldi_funs.R')
source('scripts/worland/utils.R')

# load sites and hucs
sites <- read_feather("data/decade1950plus_site_list.feather")
huc12s <- read_feather("data/huc12_list_comid.feather")

# national inventory of damns ----
item_list <- read_csv("data/basinchars/nhd_sb/basin_sb_items.csv") %>%
  filter(item=="5910de31e4b0e541a03ac983") 

houses <- sw_sb_extract(item_list$item)

# create gage housing density data frame
house_gage_df <- as.data.frame(houses$gages) %>%
  select(-matches("CAT|ACC|NODATA")) %>%
  rename(comid=COMID) %>%
  left_join(sites, by="comid") %>%
  distinct(site_no,.keep_all = T) %>%
  select(comid,site_no,everything()) %>%
  gather(variable,value,-comid,-site_no) %>%
  mutate(decade=parse_number(variable),
         decade=recode(decade,
                       '40'=1940,
                       '50'=1950,
                       '60'=1960,
                       '70'=1970,
                       '80'=1980,
                       '90'=1990,
                       '0'=2000,
                       '10'=2010),
         variable=gsub("\\d", "", variable)) %>%
  spread(variable,value) %>%
  filter(decade > 1940) 

write_feather(house_gage_df,"data/basinchars/nhd_sb/housedens_gage_df.feather")

# create huc housing density data frame
house_huc12_df <- as.data.frame(houses$hucs) %>%
  select(-matches("CAT|ACC|NODATA")) %>%
  rename(comid=COMID) %>%
  left_join(huc12s, by="comid") %>%
  distinct(huc12,.keep_all = T) %>%
  select(comid,huc12,everything()) %>%
  gather(variable,value,-comid,-huc12) %>%
  mutate(decade=parse_number(variable),
         decade=recode(decade,
                       '40'=1940,
                       '50'=1950,
                       '60'=1960,
                       '70'=1970,
                       '80'=1980,
                       '90'=1990,
                       '0'=2000,
                       '10'=2010),
         variable=gsub("\\d", "", variable)) %>%
  spread(variable,value) %>%
  filter(decade > 1940) 

write_feather(house_huc12_df,"data/basinchars/nhd_sb/housedens_huc12_df.feather")
