
library(pacman)
pacman::p_load(tidyverse,stringr,lubridate,feather,sf,proxy,lmomco)
source("scripts/worland/utils.R")

d <- read_feather("data/gage/all_gage_data.feather")

# Load daily streamflow
load("data/dvs/DV.RData")
dv_list <- as.list(DV)[d$site_no]

# randomly sample several sites
set.seed(1)
dv_sublist <- dv_list[sample(1:length(dv_list),2)]

# fdc for every gage and decade
all_direct_est <- read_feather("data/gage/all_fdc_direct_est.feather")

# CV fdc for gages and decade where observed data is available
direct_est <- read_feather("data/gage/direct_estimates.feather") %>%
  select(site_no,decade,f=variable,q_cv=nnet) %>%
  mutate(decade=as.integer(decade))

# combine and replace all_fdc with CV est
yhat <- all_direct_est %>%
  left_join(direct_est,by=c("site_no","decade","f")) %>%
  mutate(q=ifelse(!is.na(q_cv),q_cv,q),
         f=f/100) %>%
  select(site_no,decade,f,q)

test_sites <- unique(yhat$site_no)
decades <- rev(unique(yhat$decade))

all_x <- read_feather("data/gage/all_gage_covariates.feather") %>%
  mutate(decade=as.character(decade))

ref_x <- read_feather("data/gage/all_gage_data.feather")  %>%
  select(comid,site_no,decade) %>%
  left_join(all_x,by=c("comid","site_no","decade"))

est_dvlist <- list()
#for(i in 1:length(test_sites)){
for(i in 10:15){
  
  # select site
  test_site <- test_sites[i]
  
  # find reference sites for each decade
  ref_sites <- sw_find_ref(all_x,ref_x,site=test_site,distance=100)$refs
  
  Qest_all <- NULL
  for(j in 1:length(decades)) {
    
    # select decade
    test_decade <- decades[j]
    
    # EP time series for reference site
    ref_site <- ref_sites %>%
      filter(decade==test_decade) %>%
      select(ref_site) %>%
      as.character()
    
    donor_ep <- dv_list[ref_site][[1]] %>%
      filter(decade==test_decade) %>%
      mutate(ep = sw_efdc(Flow)) %>%
      select(date=Date,ep)
    
    # use estimated FDC to calculate daily flows
    est_fdc <- yhat %>%
      filter(site_no==test_site & decade==test_decade)
    
    fit <- loess(q~f,data=est_fdc,span=0.2)
    
    Q_est <- data.frame(date=donor_ep$date,
                       Q_est=round(predict(fit,donor_ep$ep),0)) 
    
    Qest_all <- rbind(Qest_all,Q_est)
  }
   
  Q_obs <- dv_list[test_site][[1]] %>%
    select(date=Date,Q_obs=Flow)
                    
  Q_all <- Qest_all %>%
    arrange(date) %>%
    left_join(Q_obs,by="date") %>%
    mutate(decade=year(floor_date(date, years(10))))
  
  est_dvlist[[i]] <- Q_all
  names(est_dvlist)[i] <- eval(test_site)
  
}

ggplot(est_dvlist[[14]]) +
  geom_line(aes(date,Q_obs)) +
  geom_line(aes(date,Q_est),color="dodgerblue",alpha=0.7) +
  facet_wrap(~decade,scales="free",ncol=1) +
  theme_bw() +
  ggtitle('Predicted streamflow',subtitle="blue=estimated, back=observed") +
  labs(y="Q")



