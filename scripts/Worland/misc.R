
library(tidyverse)
library(lubridate)
library(feather)
library(geoknife)

d <- read_feather("data/basinchars/old_sb/BASIN_CHAR_TOT.feather") 
sites <- read_csv("data/decade1950plus_site_list.csv")

por <- d %>%
  select(siteno,min_date,max_date,lat=dec_lat_va,lon=dec_long_v) %>%
  mutate(siteno = paste0("0",siteno),
         start = year(min_date),
         end = year(max_date),
         length = end-start,
         kept = ifelse(siteno %in% sites$site_no,"retained","removed")) %>%
  arrange(start,length) %>%
  na.omit() %>%
  mutate(id=rev(1:nrow(.)))

ggplot(por) + 
  geom_segment(aes(x = start, y = id, xend = end, yend = id,color=kept),alpha=0.3) +
  scale_color_manual(values=c("grey","slateblue3"),name="period of record") +
  labs(x="years",y="Individual gages") +
  geom_vline(xintercept=1950,linetype="dashed") +
  theme_bw(base_size=18) +
  theme(legend.position = c(.3, .15),
        axis.text.y=element_blank(),
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

mx <- data.frame(sites)
for (i in 1:length(dv_list)){
  date <- dv_list[[i]]$Date
  mx$max_date[i] <- as.character(max(date))
}

mx <- filter(mx,sites %in% d$site_no)
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

# plot AEP
library(lmomco)
L1 = 10^(2.13)
T2 = 0.65
T3 = 0.56
T4 = 0.36
lmr <- vec2lmom(c(L1,T2,T3,0.5*T4), lscale=FALSE)
par <- paraep4(lmr, snap.tau4=TRUE) 
FF <- seq(0,1,0.005)
Q <- qlmomco(FF, par)


FF <- c(0.00001, 0.0001, 0.001, seq(0.01, 0.99, by=0.01),
       0.999, 0.9999, 0.99999)

par1 <- list(para=c(0,1,1,1), type="aep4")
par2 <- list(para=c(0,3,1,1), type="aep4")
par3 <- list(para=c(0,5,1,1), type="aep4")

aep <- data.frame(FF) %>%
  mutate("0,1,1,1" = quaaep4(FF, par1),
         "0,3,1,1" = quaaep4(FF, par2),
         "0,5,1,1" = quaaep4(FF, par3)) %>%
  gather(parameters,value,-FF)

ggplot(aep) + 
  geom_line(aes(FF,value,linetype=parameters)) +
  scale_linetype_manual(values=c("solid","dotted","longdash")) +
  coord_cartesian(ylim=c(-10,10)) +
  xlab("Nonexceedance Probability") +
  theme_bw() +
  theme(legend.position = c(0.85,0.25), 
        legend.box.background = element_rect())

# plot activation functions
x <- seq(-10,10,0.01)

sigmoid = function(x){1/(1+exp(-x))}
tanh = function(x){2*(sigmoid(1*x))-1}
relu = function(x){sapply(x, function(z) max(0,z))}

d <- data.frame(x=x,
                sigmoid=sigmoid(x),
                tanh=tanh(x),
                relu=relu(x)) %>%
  gather(activation, value, -x) 

ggplot(d) + 
  geom_line(aes(x,value,linetype=activation)) +
  scale_linetype_manual(values=c("solid","dotted","longdash")) +
  coord_cartesian(ylim=c(-1,1),xlim=c(-5,5)) +
  labs(y="f(x)") +
  theme_bw()

d <- data.frame(x=x,
                "bias1"=sigmoid(x-2),
                "bias2"=sigmoid(x+0),
                "bias3"=sigmoid(x+2)) %>%
  gather(bias, value, -x) %>%
  mutate(bias = ifelse(bias=="bias1","sigmoid(x-2)",bias),
         bias = ifelse(bias=="bias2","sigmoid(x+0)",bias),
         bias = ifelse(bias=="bias3","sigmoid(x+2)",bias))

ggplot(d) + 
  geom_line(aes(x,value,linetype=bias)) +
  scale_linetype_manual(values=c("solid","dotted","longdash")) +
  labs(linetype="varying bias",y="f(x)") +
  theme_bw()

# possible physics based L1
d <- read_feather("data/gage/all_gage_data.feather") 

l1_area <- d %>%
  select(site_no,L1,area=tot_basin_area) %>%
  mutate(L1 = round(log10(L1)),
         area = round(log10(area))) %>%
  arrange(L1,area)

ggplot(l1_area) + geom_point(aes(L1,tot_basin_area),alpha=0.4)


