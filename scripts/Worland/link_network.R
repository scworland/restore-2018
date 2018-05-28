

library(foreign)
library(sf)
library(tidyr)
library(ggnetwork)

updated_comids <- read.dbf("data/gage/Gages_Order_PathIDs.dbf", as.is = FALSE) %>%
  select(comid=FEATUREID,site_no,huc12=HUC_12) %>%
  mutate(comid = as.character(comid)) %>%
  select(-huc12)

write_feather(updated_comids,"data/gage_list_comid.feather")

write_feather(updated_comids, "data/gage_list_comid.feather")
nhd <- st_read('data/shapefiles/paths/restore_gages_order_paths.shp',stringsAsFactors = F)

ggplot(nhd) +
  geom_sf(aes(color=totGAGE_OR),show.legend="point") +
  scale_color_viridis_c() +
  theme_bw()

site_info <- st_read('data/shapefiles/gages/restore_gages.shp',stringsAsFactors = F)

d <- read_feather("data/gage/all_gage_data.feather") %>%
  left_join(site_info,by="site_no") %>%
  select(site_no,lon,lat,L1,tot_basin_area,ppt=ppt_mean) %>%
  group_by(site_no) %>%
  summarize(area=mean(tot_basin_area),
            L1=mean(L1),
            ppt=mean(ppt)*0.00328084,
            Qn=L1/area,
            Qp=L1/ppt,
            lon=mean(lon),
            lat=mean(lat)) 

qcheck <- d %>%
  left_join(nhd,by=c("lon","lat")) %>%
  mutate(#TerminalPa = substr(TerminalPa,1,4),
    TerminalPa = as.character(TerminalPa),
    StreamOrde = as.character(StreamOrde)) %>%
  group_by(TerminalPa) %>%
  mutate(cnt = n()) %>%
  ungroup() %>%
  filter(cnt>5)

ggplot(qcheck,aes(x=totGAGE_OR,y=log10(L1),color=TerminalPa)) +
  #geom_jitter(shape=21,alpha=1,width = 0.25,height=0) +
  #geom_text(aes(label=site_no)) +
  geom_point(shape=21) +
  facet_wrap(~TerminalPa,scales="free_x") +
  labs(color="region",x="number of gages upstream from given gage",y="log10(mean streamflow)") +
  theme_bw() +
  theme(legend.position = 'none')

