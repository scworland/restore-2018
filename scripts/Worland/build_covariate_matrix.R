
library(tidyverse)
library(lubridate)
library(feather)

d <- read_feather("data/basinchars/BASIN_CHAR_TOT.feather") 

por <- d %>%
  select(siteno,min_date,max_date,lat=dec_lat_va,lon=dec_long_v) %>%
  mutate(start = year(min_date),
         end = year(max_date),
         full = ifelse(start<=1990 & end >= 2015,"1980-2010 (28%)","partial (72%)")) %>%
  arrange(full) %>%
  na.omit() %>%
  mutate(id=1:nrow(.)) 

ggplot(por) + 
  geom_segment(aes(x = start, y = id, xend = end, yend = id,color=full),alpha=0.2) +
  scale_color_manual(values=c("dodgerblue","black"),name="period of record") +
  geom_vline(xintercept=1990, linetype="dashed",color="black",size=1) +
  geom_vline(xintercept=2015, linetype="dashed",color="black",size=1) +
  labs(x="years",y="Individual gages") +
  ggtitle("Period of record for 1316 streamgages",
          subtitle="dashed lines = 1980-2010 coverage of daymet climate data") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black")) 

state <- map_data('state') 

ggplot() + geom_polygon(data=state,aes(long,lat,group=group),color = "black", fill= "white",size=1) +
  coord_fixed(1.3) + theme_bw() + xlab("") + ylab("") + 
  geom_point(data=por,aes(lon,lat,fill=full),shape=21,alpha=0.8) +
  scale_fill_manual(values=c("dodgerblue","grey95"),name="period of record") +
  theme_void() + 
  coord_fixed(1.3,xlim=c(-102,-80),ylim=c(25,38)) +
  theme(legend.position = "bottom")

# sliding window

num <- numeric()
period <- character()
for (i in 1:16) {
  start_y <- 1979+i
  end_y <- 1999+i
  
  hold <- d %>%
    select(siteno,min_date,max_date,lat=dec_lat_va,lon=dec_long_v) %>%
    mutate(start = year(min_date),
           end = year(max_date)) %>%
    filter(start<=start_y & end>= end_y)
  
  num[i] <- nrow(hold)
  period[i] <- paste0(start_y,"-",end_y)
}

por_window <- data.frame(window=period,
                         number=num) %>%
  mutate(proportion=paste0(round(number/1316,2)*100,"%"))




                 
                 
                 
                 