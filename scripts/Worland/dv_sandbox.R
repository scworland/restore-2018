
library(tidyverse)
library(lmomco)
library(lubridate)
library(feather)
source("scripts/worland/utils.R")

# Load daily streamflow
load("data/dvs/DV.RData")
sites <- ls(DV)
dv_list <- as.list(DV)
dv_list <- remove_aux(dv_list) 

# randomly sample several sites
set.seed(1)
dv_sublist <- dv_list[sample(1:length(dv_list),12)]

# create data frame
dvs <- do.call("rbind", dv_sublist) %>%
  rename(siteno=site_no,date=Date,Q=Flow,
         cd=Flow_cd)

n = length(unique(dvs$siteno))

# lmom3 <- dvs %>%
#   select(siteno,date,Q) %>%
#   mutate(Q=log(Q+0.001)) %>%
#   group_by(siteno) %>% 
#   do(lms=lmoms(.$Q, nmom=3)$lambdas) %>%
#   ungroup() %>%
#   unnest(lms) %>%
#   mutate(names = rep(c("lm1","lm2","lm3"),n)) %>%
#   spread(names,lms) %>%
#   left_join(gs,by="siteno")

fdcs <- dvs %>%
  select(siteno,date,decade,Q) %>%
  mutate(decade=as.character(decade)) %>%
  group_by(siteno,decade) %>%    
  mutate(p = sw_efdc(Q+0.001))

ggplot(fdcs) + 
  geom_line(aes(p,Q,color=decade)) + 
  facet_wrap(~siteno, scales="free_y", ncol=3) +
  scale_y_continuous(trans='log2') +
  scale_color_brewer(palette = "YlOrRd") +
  theme_dark()

# ggplot(fdcs) + 
#   geom_line(aes(p,Q)) + 
#   geom_line(aes(p2,Q),linetype="dashed",color="dodgerblue") +
#   facet_wrap(~siteno, scales="free_y") +
#   scale_y_continuous(trans='log2') +
#   theme_bw()


# prob x < value 
prob <- pnorm(q=c(28,100), mean=50, sd=20)

# quantile associated with prob
quant <- qnorm(p=0.8,mean=50,sd=20)

x <- 5
y <- 3

alpha <- rnorm(10,2,0.5)
beta <- (y-alpha)/x

ggplot() + 
  geom_point(aes(x,y)) +
  coord_cartesian(xlim=c(0,10),ylim=c(0,10)) +
  theme_bw() +
  geom_abline(intercept=alpha,slope=beta)

hold <- read_feather("data/gage/fdc_lmr_pplo.feather")
