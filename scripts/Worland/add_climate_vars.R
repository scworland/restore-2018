
library(pacman)
pacman::p_load(tidyverse,stringr,lubridate,feather,geoknife,sf,sp,dataRetrieval)

# uncomment to run from command line
# setwd("/Users/scworlan/Documents/Restore")

# source scripts
source('scripts/worland/nldi_funs.R')
source('scripts/worland/utils.R')

# load sites and hucs
sites <- read_feather("data/decade1950plus_site_list.feather")
huc12s <- read_feather("data/huc12_list_comid.feather")

# precip ----
item_list <- read_csv("data/basinchars/nhd_sb/basin_sb_items.csv") %>%
  filter(item=="57bf5c07e4b0f2f0ceb75b1b") 

precip <- sw_sb_extract(item_list$item)

# create gage precip data frame
precip_gage_df <- as.data.frame(precip$gages) %>%
  select(-matches("CAT|ACC")) %>%
  rename(comid=COMID) %>%
  left_join(sites, by="comid") %>%
  distinct(site_no,.keep_all = T) %>%
  gather(variable,value,-comid,-site_no) %>%
  mutate(year=parse_number(variable),
         decade=as.numeric(sub("\\d$", "", year))*10) %>%
  group_by(comid,site_no,decade) %>%
  summarize(ppt_mean = mean(value),
            ppt_sd = sd(value)) %>%
  ungroup() %>%
  filter(decade != 1940) %>%
  select(comid, site_no, everything())

# create huc12 precip data frame
precip_huc12_df <- as.data.frame(precip$hucs) %>%
  select(-matches("CAT|ACC")) %>%
  rename(comid=COMID) %>%
  left_join(huc12s, by="comid") %>%
  distinct(huc12,.keep_all = T) %>%
  gather(variable,value,-comid,-huc12) %>%
  mutate(year=parse_number(variable),
         decade=as.numeric(sub("\\d$", "", year))*10) %>%
  group_by(comid,huc12,decade) %>%
  summarize(ppt_mean = mean(value),
            ppt_sd = sd(value)) %>%
  ungroup() %>%
  filter(decade != 1940) %>%
  select(comid, huc12, everything())

# temperature ----
item_list <- read_csv("data/basinchars/nhd_sb/basin_sb_items.csv") %>%
  filter(item=="5787ea72e4b0d27deb377b6d") 

temp <- sw_sb_extract(item_list$item)

# create gage temp data frame
temp_gage_df <- as.data.frame(temp$gages) %>%
  select(-matches("CAT|ACC")) %>%
  rename(comid=COMID) %>%
  left_join(sites, by="comid") %>%
  distinct(site_no,.keep_all = T) %>%
  gather(variable,value,-comid,-site_no) %>%
  mutate(year=parse_number(variable),
         decade=as.numeric(sub("\\d$", "", year))*10) %>%
  group_by(comid,site_no,decade) %>%
  summarize(temp_mean = mean(value),
            temp_sd = sd(value)) %>%
  ungroup() %>%
  filter(decade != 1940) %>%
  select(comid, site_no, everything())

write_feather(temp_gage_df,"data/basinchars/nhd_sb/temp_gage_df.feather")

# create huc12 temp data frame
temp_huc12_df <- as.data.frame(temp$hucs) %>%
  select(-matches("CAT|ACC")) %>%
  rename(comid=COMID) %>%
  left_join(huc12s, by="comid") %>%
  distinct(huc12,.keep_all = T) %>%
  gather(variable,value,-comid,-huc12) %>%
  mutate(year=parse_number(variable),
         decade=as.numeric(sub("\\d$", "", year))*10) %>%
  group_by(comid,huc12,decade) %>%
  summarize(temp_mean = mean(value),
            temp_sd = sd(value)) %>%
  ungroup() %>%
  filter(decade != 1940) %>%
  select(comid, huc12, everything())

write_feather(temp_huc12_df,"data/basinchars/nhd_sb/temp_huc12_df.feather")

# combine climate data ----
# precip_gage_df <- read_feather("data/basinchars/nhd_sb/precip_gage_df.feather")
# precip_huc12_df <- read_feather("data/basinchars/nhd_sb/precip_huc12_df.feather")
# temp_gage_df <- read_feather("data/basinchars/nhd_sb/temp_gage_df.feather")
# temp_huc12_df <- read_feather("data/basinchars/nhd_sb/temp_huc12_df.feather")

climate_gage_df <- precip_gage_df %>%
  select(comid,decade,ppt_mean,ppt_sd) %>%
  left_join(temp_gage_df,by=c("comid","decade")) %>%
  distinct(site_no,decade,.keep_all = T) %>%
  select(comid, site_no, everything())

write_feather(climate_gage_df,"data/basinchars/nhd_sb/climate_gage_df.feather")

climate_huc12_df <- precip_huc12_df %>%
  select(comid,decade,ppt_mean,ppt_sd) %>%
  left_join(temp_huc12_df,by=c("comid","decade")) %>%
  distinct(huc12,decade,.keep_all = T) %>%
  select(comid, huc12, everything())

write_feather(climate_huc12_df,"data/basinchars/nhd_sb/climate_huc12_df.feather")
