
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


# Kmeans for coordinates
site_info <- st_read('data/shapefiles/gages/restore_gages.shp') %>%
  mutate(site_no=as.character(site_no))

km_coords <- select(site_info,lat,lon) %>%
  st_set_geometry(NULL) %>%
  do(data.frame(., kclust = kmeans(as.matrix(.),centers=6)$cluster)) %>%
  mutate(site_no = site_info$site_no,
         kclust = as.character(kclust)) %>%
  select(-lat,-lon) %>%
  left_join(site_info, by="site_no")

ggplot(km_coords) + 
  geom_sf(aes(shape=kclust),show.legend="point") + 
  labs(x="Longitude",y="Latitude") +
  #coord_sf(datum = NA) +
  theme_bw() 

# function to find reference site
site_info <- st_read('data/shapefiles/gages/restore_gages.shp',stringsAsFactors = F) 

d <- read_feather("data/gage/all_gage_data.feather") 

all_x <- read_feather("data/gage/all_gage_covariates.feather") %>%
  mutate(decade=as.character(decade))

ref_x <- d %>%
  select(comid,site_no,decade) %>%
  left_join(all_x,by=c("comid","site_no","decade"))

ref_list <- find_ref(all_x,ref_x,site="08023080",distance=150)

refs_plot <- ref_list$subs %>%
  gather(type,site_no,-decade) %>%
  left_join(site_info,by="site_no")

sites_plot <- site_info %>%
  filter(!site_no %in% refs_plot$site_no)

ggplot() +
  geom_sf(data=sites_plot, color="grey", alpha=0.5, size=0.8) +
  geom_sf(data=refs_plot,aes(shape=type,fill=type),show.legend="point") +
  scale_shape_manual(values=c(3,21,21)) +
  scale_fill_manual(values=c("black","orange","dodgerblue")) +
  facet_wrap(~decade) +
  theme_bw() +
  theme(legend.position = 'top')




