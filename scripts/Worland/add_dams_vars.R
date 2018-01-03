
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
  filter(item=="58c301f2e4b0f37a93ed915a") 

dams <- sw_sb_extract(item_list$item)

# create gage dam data frame
dam_gage_df <- as.data.frame(dams$gages) %>%
  select(-matches("CAT|ACC|DENS")) %>%
  rename(comid=COMID) %>%
  left_join(sites, by="comid") %>%
  distinct(site_no,.keep_all = T) %>%
  select(comid,site_no,TOT_NDAMS1930:TOT_MAJOR2010) %>%
  gather(variable,value,-comid,-site_no) %>%
  mutate(decade=parse_number(variable),
         variable=gsub("\\d", "", variable)) %>%
  spread(variable,value) %>%
  filter(decade > 1940) 

write_feather(dam_gage_df,"data/basinchars/nhd_sb/dam_gage_df.feather")

# create create huc12 dam data frame
dam_huc12_df <- as.data.frame(dams$hucs) %>%
  select(-matches("CAT|ACC|DENS")) %>%
  rename(comid=COMID) %>%
  left_join(huc12s, by="comid") %>%
  distinct(huc12,.keep_all = T) %>%
  select(comid,huc12,TOT_NDAMS1930:TOT_MAJOR2010) %>%
  gather(variable,value,-comid,-huc12) %>%
  mutate(decade=parse_number(variable),
         variable=gsub("\\d", "", variable)) %>%
  spread(variable,value) %>%
  filter(decade > 1940) 

write_feather(dam_huc12_df,"data/basinchars/nhd_sb/dam_huc12_df.feather")
