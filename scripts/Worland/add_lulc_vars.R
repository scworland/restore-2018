library(pacman)
pacman::p_load(tidyverse,stringr,lubridate,feather,geoknife,sf,sp,dataRetrieval,corrr)

# uncomment to run from command line
setwd("/Users/scworlan/Documents/Restore")

# source scripts
source('scripts/worland/nldi_funs.R')
source('scripts/worland/utils.R')

# load sites and hucs
sites <- read_feather("data/decade1950plus_site_list.feather")
huc12s <- read_feather("data/huc12_list_comid.feather")

# historic lulc, 1940-1990
historic_class_link <- read_csv("data/basinchars/nhd_sb/historic_lulc_classes.csv") %>%
  mutate(variable=sub('.*_', '', variable))

item_list <- read_csv("data/basinchars/nhd_sb/basin_sb_items.csv") %>%
  filter(item=="58cbeef2e4b0849ce97dcd61")

historic_lulc <- sw_sb_extract(item_list$item)

historic_lulc_gages_df <- as.data.frame(historic_lulc$gages) %>%
  select(-matches("NODATA")) %>%
  rename(comid=COMID) %>%
  left_join(sites, by="comid") %>%
  distinct(site_no,.keep_all = T) %>%
  select(comid,site_no,everything()) %>%
  gather(variable,value,-comid,-site_no) %>%
  mutate(decade = sub('.*TOT_',"",variable) %>%
           sub('_.*',"",.) %>%
           parse_number(.) %>%
           paste0("19",.) %>%
           as.numeric(.),
         variable=sub('.*_', '', variable)) %>%
  left_join(historic_class_link, by="variable") %>%
  filter(lulc_class !="Intentionally left blank",
         lulc_class !="Mining") %>%
  select(-variable) %>%
  distinct(site_no,lulc_class,decade,.keep_all = T) %>%
  select(comid,site_no,decade,variable=lulc_class,value) %>%
  spread(variable,value)

historic_lulc_huc12_df <- as.data.frame(historic_lulc$hucs) %>%
  select(-matches("NODATA")) %>%
  rename(comid=COMID) %>%
  left_join(huc12s, by="comid") %>%
  distinct(huc12,.keep_all = T) %>%
  select(comid,huc12,everything()) %>%
  gather(variable,value,-comid,-huc12) %>%
  mutate(decade = sub('.*TOT_',"",variable) %>%
           sub('_.*',"",.) %>%
           parse_number(.) %>%
           paste0("19",.) %>%
           as.numeric(.),
         variable=sub('.*_', '', variable)) %>%
  left_join(historic_class_link, by="variable") %>%
  filter(lulc_class !="Intentionally left blank",
         lulc_class !="Mining") %>%
  select(-variable) %>%
  distinct(huc12,lulc_class,decade,.keep_all = T) %>%
  select(comid,huc12,decade,variable=lulc_class,value) %>%
  spread(variable,value)


# 2001 lulc
class_link <- read_csv("data/basinchars/nhd_sb/lulc_classes.csv") %>%
  mutate(variable=as.character(variable))

item_list <- read_csv("data/basinchars/nhd_sb/basin_sb_items.csv") %>%
  filter(item=="5761b67de4b04f417c2d30ae") 

lulc2001 <- sw_sb_extract(item_list$item)

lulc2001_gages_df <- as.data.frame(lulc2001$gages)  %>%
  select(-matches("NODATA")) %>%
  rename(comid=COMID) %>%
  left_join(sites, by="comid") %>%
  distinct(site_no,.keep_all = T) %>%
  select(comid,site_no,everything()) %>%
  gather(variable,value,-comid,-site_no) %>%
  mutate(decade = 2000,
         variable=sub('.*_', '', variable)) %>%
  left_join(class_link, by="variable") %>%
  select(-variable) %>%
  distinct(site_no,lulc_class,decade,.keep_all = T) %>%
  spread(lulc_class,value) %>%
  mutate(developed = developed_low,developed_medium,developed_high) %>%
  select(-c(developed_open,developed_low,developed_medium,developed_high)) 

lulc2001_huc12_df <- as.data.frame(lulc2001$hucs)  %>%
  select(-matches("NODATA")) %>%
  rename(comid=COMID) %>%
  left_join(huc12s, by="comid") %>%
  distinct(huc12,.keep_all = T) %>%
  select(comid,huc12,everything()) %>%
  gather(variable,value,-comid,-huc12) %>%
  mutate(decade = 2000,
         variable=sub('.*_', '', variable)) %>%
  left_join(class_link, by="variable") %>%
  select(-variable) %>%
  distinct(huc12,lulc_class,decade,.keep_all = T) %>%
  spread(lulc_class,value) %>%
  mutate(developed = developed_low,developed_medium,developed_high) %>%
  select(-c(developed_open,developed_low,developed_medium,developed_high)) 

# 2011 lulc
item_list <- read_csv("data/basinchars/nhd_sb/basin_sb_items.csv") %>%
  filter(item=="5761bad4e4b04f417c2d30c5")

lulc2011 <- sw_sb_extract(item_list$item,type="ACC")

lulc2011_gages_df <- as.data.frame(lulc2011$gages)  %>%
  select(-matches("NODATA")) %>%
  rename(comid=COMID) %>%
  left_join(sites, by="comid") %>%
  distinct(site_no,.keep_all = T) %>%
  select(comid,site_no,everything()) %>%
  gather(variable,value,-comid,-site_no) %>%
  mutate(decade = 2010,
         variable=sub('.*_', '', variable)) %>%
  left_join(class_link, by="variable") %>%
  select(-variable) %>%
  distinct(site_no,lulc_class,decade,.keep_all = T) %>%
  spread(lulc_class,value) %>%
  mutate(developed = developed_low,developed_medium,developed_high) %>%
  select(-c(developed_open,developed_low,developed_medium,developed_high)) 

lulc2011_huc12_df <- as.data.frame(lulc2011$hucs)  %>%
  select(-matches("NODATA")) %>%
  rename(comid=COMID) %>%
  left_join(huc12s, by="comid") %>%
  distinct(huc12,.keep_all = T) %>%
  select(comid,huc12,everything()) %>%
  gather(variable,value,-comid,-huc12) %>%
  mutate(decade = 2010,
         variable=sub('.*_', '', variable)) %>%
  left_join(class_link, by="variable") %>%
  select(-variable) %>%
  distinct(huc12,lulc_class,decade,.keep_all = T) %>%
  spread(lulc_class,value) %>%
  mutate(developed = developed_low,developed_medium,developed_high) %>%
  select(-c(developed_open,developed_low,developed_medium,developed_high)) 

# combine everything
lulc_gage_df <- bind_rows(historic_lulc_gages_df,lulc2001_gages_df,lulc2011_gages_df)
write_feather(lulc_gage_df,"data/basinchars/nhd_sb/lulc_gage_df.feather")

lulc_huc12_df <- bind_rows(historic_lulc_huc12_df,lulc2001_huc12_df,lulc2011_huc12_df) %>%
  filter(decade != 1940)

write_feather(lulc_huc12_df,"data/basinchars/nhd_sb/lulc_huc12_df.feather")

