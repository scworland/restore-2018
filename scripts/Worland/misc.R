
library(tidyverse)
library(lubridate)
library(feather)
library(geoknife)

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


# OLD DV sandbox

# load R environment
load("data/dvs/DV.RData")

# extract site numbers
sites <- ls(DV)

# load list of dvs
dv_list <- as.list(DV)

# LULC files
pth = "data/basinchars/sciencebase"
lulc_files <- list.files(path=pth, pattern = glob2rx("*LULC*.feather"),full.names = T)
temp_files <- list.files(path=pth, pattern = glob2rx("*TAV*.feather"),full.names = T)[2:13]
ppt_files <- list.files(path=pth, pattern = glob2rx("*PPT*.feather"),full.names = T)

temp <- read_months(temp_files,"temp") %>%
  reduce(left_join, by = "siteno")

ppt <- read_months(ppt_files,"ppt") %>%
  reduce(left_join, by = "siteno")

lulc <- read_years(lulc_files) %>%
  reduce(left_join, by = "siteno") %>%
  select(matches('Cult|Dev|Dec|Ever|Hay|Mixed')) %>%
  select(contains('2004'))

# load basin chars
bfi <- read_feather("data/basinchars/BFI.feather")
et <- read_feather("data/basinchars/ET.feather")
runoff <- read_feather("data/basinchars/RUN7100.feather")
tav <- read_feather("data/basinchars/TAV7100_ANN.feather")

df_list <- list(bfi,et,runoff,tav) 

X <- df_list %>%
  reduce(left_join, by = "siteno") %>%
  select(siteno,contains("CAT")) %>%
  left_join(temp) %>% 
  left_join(ppt) %>%
  left_join(lulc)

# PCA regression
d <- read_feather("data/gage/lmom3.feather")

X <- select(d,siteno,cat_area_sqkm:tot_wildfire_2011) %>%
  gather(variable,value,-siteno) %>%
  group_by(siteno) %>%
  filter(value > 0) %>% 
  ungroup() %>%
  spread(variable, value) %>%
  select_if(~!any(is.na(.))) %>%
  select(-siteno)


pcs <- prcomp(X, scale=T)

# screeplot(pcs,npcs=20,type="lines")

Xpc <- cbind(select(d,siteno:comid),pcs$x[,1:8])

lm1 <- lm(lm1 ~., data=select(Xpc,lm1,PC1:PC8))
lm2 <- lm(lm2 ~., data=select(Xpc,lm2,PC1:PC8))
lm3 <- lm(lm3 ~., data=select(Xpc,lm3,PC1:PC8))

resid <- data.frame(lm1r = residuals(lm1),
                    lm2r = residuals(lm2),
                    lm3r = residuals(lm3))


# download FTP data
library(curl)
url = "ftp://ftpext.usgs.gov/pub/er/wi/middleton/dblodgett/catchment_chars/rds/"
h = new_handle(dirlistonly=TRUE)
con = curl(url, "r", h)
tbl = read.table(con, stringsAsFactors=TRUE, fill=TRUE)
close(con)
head(tbl)

urls <- paste0(url, tbl[1,1])
fls = basename(urls)
curl_fetch_disk(urls[1], fls[1])

library(readr)
library(sbtools)
hold <- read_delim("data/basinchars/nhd_sb/SOHL60_TOT_CONUS.txt",",")

sw_get_files <- function(token){
  
  files <- item_list_files(token) %>%
    filter(grepl("TOT",fname))
  
  for (i in 1:nrow(files)){
    pth <- file.path("data/basinchars/nhd_sb",files$fname[i])
    
    item_file_download(token, 
                       names = files$fname[i], 
                       destinations = pth)
    
    unzip(pth,exdir="data/basinchars/nhd_sb")
  }
}

sw_get_files('58cbeef2e4b0849ce97dcd61')


