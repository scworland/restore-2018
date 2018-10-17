
sw_predict_fdcs <- function(final_fdc_model,gage_covariates,huc12_covariates,mnn_quantile_est) {
  
  model <- final_fdc_model
  
  # combine huc12 and all gage covariates
  all_covariates <- gage_covariates %>%
    bind_rows(huc12_covariates) %>% # bind rows
    select(comid,site_no,decade,everything()) %>%
    distinct(comid,decade,.keep_all=T) # drop duplicate comids
  
  # prepare covariates for model predictions
  all_X <- all_covariates %>%
    select(major:flood_storage, -basin_area) %>%
    mutate_all(funs(as.numeric(as.factor(.)))) %>%
    mutate_all(funs(as.vector(scale(.)))) %>%
    select_if(~!any(is.na(.))) %>%
    as.matrix()
  
  cv_yhat <- mnn_quantile_est$est_obs %>%
    select(site_no,decade,nep=variable,cv_q=nnet, obs_q=obs)
  
  # predict FDCs for every site
  predictions <- predict(model, all_X)
  
  # check number of violations
  viol_all <- predictions %>%
    data.frame() %>%
    setNames(colnames(Y)) %>%
    mutate(comid=all_covariates$comid, 
           decade=as.character(all_covariates$decade),
           area = all_covariates$basin_area) %>%
    gather(f,q,-comid,-decade,-area) %>%
    mutate(q=10^(q)*area) %>%
    group_by(comid,decade) %>%
    summarize(viol = sum(cummax(q) != q))
  
  viol <- sum(yhat$viol)/(nrow(yhat) * 27)

  yhat <- predictions %>%
    data.frame() %>%
    setNames(colnames(Y)) %>%
    mutate(comid=all_covariates$comid, # add comid
           site_no=all_covariates$site_no, # add site_no
           huc12=all_covariates$huc12, # add huc12s
           decade=as.character(all_covariates$decade), # add decade
           lon=all_covariates$dec_long_va, # add lon
           lat=all_covariates$dec_lat_va,
           area=all_covariates$basin_area) %>% # add lat
    gather(nep,q,-comid,-site_no,-huc12,-decade,-lon,-lat,-area) %>% # convert wide --> long
    mutate(nep = as.numeric(substring(nep, 2))) %>%
    mutate(q = round(10^(q)*area,2)) %>%
    left_join(cv_yhat, by = c("site_no","decade","nep")) %>%
    mutate(site_no = ifelse(is.na(cv_q), NA, site_no)) %>%
    arrange(comid,decade)
  
  return(yhat)
}

# write_feather(yhat,"data/gage/all_fdc_direct_est.feather")
