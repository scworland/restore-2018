
library(pacman)
pacman::p_load(tidyverse,stringr,lubridate,feather,sf,proxy,lmomco)
source("scripts/worland/utils.R")

d <- read_feather("data/gage/all_gage_data.feather") %>%
  filter(decade != 2010)

# Load daily streamflow
load("data/dvs/DV.RData")
# sites <- ls(DV)
dv_list <- as.list(DV)[d$site_no]

# randomly sample several sites
set.seed(1)
dv_sublist <- dv_list[sample(1:length(dv_list),2)]

# create data frame
dvs <- do.call("rbind", dv_sublist) %>%
  rename(date=Date,Q=Flow,cd=Flow_cd)

n = length(unique(dvs$siteno))

estimate_dv <- function(lmoms,basinchars,dv_list){
  for (i in 1:nrow(lmoms)) {
    
    print(i)
    
    # extract site and decade
    test_site <- lmoms$site_no[i]
    test_decade <- lmoms$decade[i]
    
    # subset covariates for decade
    X <- select(basinchars,site_no,decade,ppt_mean:tot_rdx) %>%
      mutate_at(vars(-site_no,-decade),funs(as.vector(scale(.)))) %>%
      select_if(~!any(is.na(.))) %>%
      filter(decade==test_decade)
    
    # find donor site
    dists <- proxy::dist(X[-i,-c(1,2)],X[i,-c(1,2)])
    donor_site <- X$site_no[which(dists==min(dists))+1]
    
    # generate EPs from donor site
    donor_ep <- dv_list[donor_site][[1]] %>%
      filter(decade==test_decade) %>%
      mutate(ep = sw_efdc(Flow)) %>%
      select(date=Date,ep)
    
    # parameterize distribution using lmoments
    test_lmoms <- lmoms %>%
      filter(site_no==test_site & decade==test_decade) %>%
      select(-site_no,-decade)
    
    test_par <- lmom2par(vec2lmom(as.numeric(test_lmoms)),type="pe3")
    test_Q <- qlmomco(donor_ep$ep, test_par)
    
    est_obs[[i]] <- dv_list[test_site][[1]] %>%
      filter(decade==test_decade) %>%
      select(site_no,date=Date,obs_Q=Flow) %>%
      mutate(est_Q=test_Q) %>%
      gather(variable,value,-site_no,-date)
  }

result <- bind_rows(est_obs)
return(result)

}

# example
estdvs <- estimate_dv(pred_lmoms[1:15,],d,dv_list) 

nse <- estdvs %>%
  spread(variable,value) %>%
  group_by(site_no) %>%
  summarize(nse = NSE(est_Q,obs_Q))

ggplot(estdvs) + 
  geom_line(aes(date,value,color=variable)) +
  scale_color_manual(values=c("darkgrey","dodgerblue")) +
  facet_wrap(~site_no) +
  theme_bw()








