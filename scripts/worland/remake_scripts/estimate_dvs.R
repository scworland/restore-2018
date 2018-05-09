
sw_estimate_dv <- function(dv_list,gage_all,all_est_fdc,gage_covariates,huc12_covariates){
  
  # load all estimated FDCs
  yhat <- all_est_fdc %>%
    mutate(f=f/100) 
  
  # prepare gage_covariates to merge
  gage_x <- gage_covariates %>%
    mutate(decade=as.character(decade),
           huc12=as.character(huc12)) %>%
    select(comid,site_no,huc12,lon,lat,decade,everything())
  
  # prepare huc12_covariates to merge 
  huc12_x <- huc12_covariates %>%
    mutate(decade=as.character(decade),
           huc12=as.character(huc12)) %>%
    select(comid,huc12,lon,lat,decade,everything()) 
  
  # combine gage and huc x and keep distinct comid+decade
  all_x <- bind_rows(gage_x,huc12_x) %>%
    distinct(comid,decade,.keep_all=T)
  
  # basin chars for only index basins
  index_x <- gage_all %>%
    select(comid,site_no,decade) %>%
    left_join(gage_x,by=c("comid","site_no","decade"))
  
  test_sites <- unique(yhat$comid)
  decades <- rev(unique(yhat$decade))
  est_dvlist <- list()
  #for(i in 1:length(test_sites)){
  for(i in 1:10){
    
    print(paste0("Predicting dvs for ",i," out of ", length(test_sites), " locations"))
    
    # select site
    test_site <- test_sites[i]
    
    # find reference sites for each decade
    index_sites <- sw_find_index(all_x,index_x,site=test_site,distance=100)$index_sites
    
    Qest_all <- NULL
    for(j in 1:length(decades)) {
      
      # select decade
      test_decade <- decades[j]
      
      # Find site_no associated with comid
      index_site <- index_sites %>%
        filter(decade==test_decade) %>%
        rename(comid=index_site) %>%
        left_join(select(gage_all,comid,site_no,decade),
                  by=c("comid","decade")) %>%
        select(site_no) %>%
        as.character()
      
      donor_ep <- dv_list[index_site][[1]] %>%
        filter(decade==test_decade) %>%
        mutate(ep = sw_efdc(Flow)) %>%
        select(date=Date,ep)
      
      # use estimated FDC to calculate daily flows
      est_fdc <- yhat %>%
        filter(comid==test_site & decade==test_decade)
      
      fit <- loess(q~f,data=est_fdc,span=0.2)
      
      Q_est <- data.frame(date=donor_ep$date,
                          Q_est=round(predict(fit,donor_ep$ep),0)) 
      
      Qest_all <- rbind(Qest_all,Q_est)
    }
    
    test_siteno <- unique(gage_all$site_no[test_site==gage_all$comid])
    
    if(identical(test_siteno, character(0))){
      
      Q_obs <- data.frame(date=Q_est$date, Q_obs=NA) 
      
    }else{
      
      Q_obs <- dv_list[test_siteno][[1]] %>%
        select(date=Date,Q_obs=Flow)
      
    }
    
    Q_all <- Qest_all %>%
      arrange(date) %>%
      left_join(Q_obs,by="date") %>%
      mutate(decade=year(floor_date(date, years(10))))
    
    
    est_dvlist[[i]] <- Q_all
    names(est_dvlist)[i] <- eval(test_site)
    
  }
  
  return(est_dvlist)
}

# ggplot(est_dvlist[[10]]) +
#   geom_line(aes(date,Q_obs)) +
#   geom_line(aes(date,Q_est),color="dodgerblue",alpha=0.7) +
#   facet_wrap(~decade,scales="free_x",ncol=2) +
#   theme_bw() +
#   # scale_y_log10() +
#   ggtitle('Predicted streamflow for site 02329500',
#           subtitle="blue=estimated, black=observed") +
#   labs(y="Q")
# 
# hold <- est_dvlist[[10]] %>%
#   select(Q_obs,Q_est) %>%
#   gather(variable,value) %>%
#   group_by(variable) %>%
#   summarize(mn = mean(value),
#             std = sd(value),
#             min = min(value),
#             max = max(value))
#             
# write_csv(all_x,"data/gage/gage_basin_characteristics.csv")
# 
# all_huc12 <- read_feather("data/huc12/all_huc12_covariates.feather") %>%
#   mutate(nodat = ifelse(bedperm=="nodata",1,0))
# 
# ggplot(all_huc12) + geom_point(aes(x=lon,y=lat,color=as.character(nodat)))
