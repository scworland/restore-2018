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

# create immutable basin characteristic matrix
item_list <- read_csv("data/basinchars/nhd_sb/basin_sb_items.csv") %>%
  filter(type=="immutable") 

items <- item_list$item
group <- item_list$grouping_var
sb_type <- item_list$sb_type

immutable_chars <- list()
immutable_gages <- list()
immutable_hucs <- list()
for (i in 1:length(items)){
  immutable_chars[[i]] <- sw_sb_extract(items[i],type=sb_type[i],group=group[i])
  immutable_gages[[i]] <- immutable_chars[[i]]$gages
  immutable_hucs[[i]] <- immutable_chars[[i]]$hucs
  print(paste0("completed ",item_list$description[i]," = ",i, "/13"))
}


immutable_gage_df <- as.data.frame(immutable_gages) %>%
  rename(COMID=comid) %>%
  select(-contains("comid",ignore.case = F)) %>%
  rename(comid=COMID) %>%
  left_join(sites, by="comid") %>%
  distinct(site_no, .keep_all = T) %>%
  select(-matches("cat|acc|nodata|tot_s|x1")) %>%
  select(comid, site_no, everything())

write_feather(immutable_gage_df,"data/basinchars/nhd_sb/immutable_gage.feather")

immutable_huc12_df <- as.data.frame(immutable_hucs) %>%
  rename(COMID=comid) %>%
  select(-contains("comid",ignore.case = F)) %>%
  rename(comid=COMID) %>%
  left_join(huc12s, by="comid") %>%
  distinct(huc12, .keep_all = T) %>%
  select(-matches("cat|acc|nodata|tot_s|x1")) %>%
  select(comid, huc12, everything())

write_feather(immutable_huc12_df,"data/basinchars/nhd_sb/immutable_huc12.feather")

