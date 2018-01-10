
library(pacman)
pacman::p_load(tidyverse,stringr,lubridate,feather,sf)
source("scripts/worland/utils.R")

d <- read_feather("data/gage/all_gage_data.feather")

site_info <- st_read('data/shapefiles/gages/restore_gages.shp') %>%
  mutate(site_no=as.character(site_no))

Y <- select(d,L1:T6) %>%
  as.matrix()

X <- select(d,site_no,decade,ppt_mean:tot_rdx) %>%
  mutate_at(vars(-c(site_no,decade)),funs(as.vector(scale(.)))) %>%
  select_if(~!any(is.na(.))) %>%
  arrange(decade)

km <- X %>%
  select(-site_no) %>%
  group_by(decade) %>%
  do(data.frame(., kclust = kmeans(as.matrix(.),centers=5)$cluster)) %>%
  ungroup() %>%
  mutate(site_no = X$site_no,
         kclust = as.character(kclust)) %>%
  select(site_no,decade,kclust) %>%
  left_join(site_info, by="site_no")

ggplot(km) + 
  geom_sf(aes(color=kclust)) + 
  facet_wrap(~decade, ncol=2) + 
  theme_minimal() +
  scale_color_viridis_d()

site = "02294650"
km_group <- obs_in_cluster(site,km) %>% 
  left_join(site_info, by="site_no")

ggplot() + 
  geom_sf(data=km,color="grey",alpha=0.5) + 
  geom_sf(data=km_group,color="dodgerblue",alpha=0.7) +
  geom_sf(data=filter(km,site_no==site),shape=23,fill="red") +
  facet_wrap(~decade, ncol=3) + 
  ggtitle("Decadal K-clusters for site #02294650",
          subtitle="Site indicated by red diamond. Based on decadal Kmeans with 5 centers calculated from 34 basin characteristics") +
  theme_bw() 

