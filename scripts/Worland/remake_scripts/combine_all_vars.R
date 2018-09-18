
sw_combine_gage <- function(sites,immutable_vars,dams,housing_density,climate,lulc){
  
  filter <- dplyr::filter
  select <- dplyr::select
  
  ## variables to be combined
  immutable_gage_df <- immutable_vars$immutable_gage_df 
  house_gage_df <- housing_density$house_gage_df
  dam_gage_df <- dams$dam_gage_df 
  climate_gage_df <- climate$climate_gage_df
  lulc_gage_df <- lulc$lulc_gage_df
  
  # used to convert acre-ft to km^2-m
  val <- 0.3048*(1/247.104393047)
  
  # list of sites impacted by edwards aquifer
  # edwards_sites <- c("08156800","08155300","08155400","08181400","08184000",
  #                    "08185000","08190500","08197500","08198500","08200700",
  #                    "08202700")
  
  # list of comids impacted by edwards aquifer
  edwards_comids <- c("5781337","5781711","5781703","10835030","7851041","7850619",
                      "7876116","10646917","10646003","10655659","10655687")
  
  # combine and filter covariates
  all_gage_covariates <- list(house_gage_df,dam_gage_df,climate_gage_df,lulc_gage_df) %>%
    reduce(left_join, by = c("comid","site_no","decade")) %>%
    arrange(comid,decade) %>%
    left_join(immutable_gage_df, by = c("comid","site_no")) %>%
    left_join(sites,by=c("comid","site_no")) %>%
    rename_all(tolower) %>%
    select(comid,site_no,huc12,lon,lat,acc_hdens:woody_wetland,
           acc_bfi:acc_twi,acc_basin_area:acc_rdx,everything()) %>%
    filter(decade %in% c(1950,1960,1970,1980,1990,2000)) %>%
    # below is WHA and SCWs covariate filter
    mutate(acc_nid_storage = acc_nid_storage*val,
           acc_norm_storage = acc_norm_storage*val,
           flood_storage = abs((acc_nid_storage - acc_norm_storage)/acc_basin_area),
           # assign 1990 flood storage to 2000 due to error
           flood_storage = ifelse(site_no == "02295420" & decade == "2000",0.12216912,flood_storage),
           acc_basin_slope = ifelse(acc_basin_slope==0,0.02,acc_basin_slope),
           soller.1 = ifelse(site_no=="02359315",soller,soller.1),
           aquifers.1 = ifelse(site_no=="02359315",aquifers,aquifers.1),
           physio.1 = ifelse(site_no=="02359315",physio,physio.1),
           aquifers = ifelse(aquifers == "cat_aq_nodata","nodata",aquifers),
           soller = ifelse(soller == "cat_soller_nodata","nodata",soller),
           aquifers.1 = ifelse(aquifers.1 == "cat_aq_nodata","nodata",aquifers.1),
           soller.1 = ifelse(soller.1 == "cat_soller_nodata","nodata",soller.1),
           ed_rch_zone = ifelse(comid %in% edwards_comids, "1", "0")) %>%
    select(-acc_elev_max,-acc_elev_min,-runoff_mean,-runoff_sd,-acc_hdens,
           -perennial_ice_snow,-acc_stream_slope) %>%
    rename(dec_lat_va=lat,dec_long_va=lon,cat_soller=soller.1,cat_physio=physio.1,
           cat_aquifers=aquifers.1,cat_ecol3=ecol3.1,length_km=acc_stream_length) %>%
    drop_na() %>%
    set_names(gsub(x = names(.), pattern = "acc_", replacement = ""))
  
    #names(all_gage_covariates) <- gsub(x = names(all_gage_covariates), pattern = "acc_", replacement = "") 
  
  return(all_gage_covariates)
}

sw_combine_huc12 <- function(huc12s,immutable_vars,dams,housing_density,climate,lulc,gage_covariates){
  
  ## variables to be combined
  immutable_huc12_df <- immutable_vars$immutable_huc12_df 
  house_huc12_df <- housing_density$house_huc12_df
  dam_huc12_df <- dams$dam_huc12_df 
  climate_huc12_df <- climate$climate_huc12_df
  lulc_huc12_df <- lulc$lulc_huc12_df
  
  # used to convert acre-ft to km^2/m
  val <- 0.3048*(1/247.104393047)
  
  # list of comids impacted by edwards aquifer
  edwards_comids <- c("5781337","5781711","5781703","10835030","7851041","7850619",
                      "7876116","10646917","10646003","10655659","10655687")
  
  all_huc12_covariates <- list(house_huc12_df,dam_huc12_df,climate_huc12_df,lulc_huc12_df) %>%
    reduce(left_join, by = c("comid","huc12","decade")) %>%
    arrange(comid,decade) %>%
    left_join(immutable_huc12_df, by = c("comid","huc12")) %>%
    left_join(huc12s, by = c("comid","huc12")) %>%
    rename_all(tolower) %>%
    select(comid,huc12,lon,lat,acc_hdens:woody_wetland,
           acc_bfi:acc_twi,acc_basin_area:acc_rdx,everything()) %>%
    filter(decade %in% c(1950,1960,1970,1980,1990,2000)) %>%
    # below is WHA and SCWs covariate filter
    mutate(acc_nid_storage = acc_nid_storage*val,
           acc_norm_storage = acc_norm_storage*val,
           flood_storage = abs((acc_nid_storage - acc_norm_storage)/acc_basin_area),
           flood_storage = ifelse(is.na(flood_storage),0,flood_storage),
           acc_basin_slope = ifelse(acc_basin_slope==0,0.02,acc_basin_slope),
           aquifers = ifelse(aquifers == "cat_aq_nodata","nodata",aquifers),
           soller = ifelse(soller == "cat_soller_nodata","nodata",soller),
           aquifers.1 = ifelse(aquifers.1 == "cat_aq_nodata","nodata",aquifers.1),
           soller.1 = ifelse(soller.1 == "cat_soller_nodata","nodata",soller.1),
           ed_rch_zone = ifelse(comid %in% edwards_comids, "1", "0")) %>%
    select(-acc_elev_max,-acc_elev_min,-runoff_mean,-runoff_sd,-acc_hdens,
           -perennial_ice_snow,-acc_stream_slope) %>%
    rename(dec_lat_va=lat,dec_long_va=lon,cat_soller=soller.1,cat_physio=physio.1,
           cat_aquifers=aquifers.1,cat_ecol3=ecol3.1,length_km=acc_stream_length) %>%
    drop_na() %>%
    set_names(gsub(x = names(.), pattern = "acc_", replacement = "")) %>%
    filter_all(all_vars(!grepl("nodata",.)))
    #mutate_at(vars(),funs(sub('.*\\_', '', .)))
  
  # names(all_huc12_covariates) <- gsub(x = names(all_huc12_covariates), pattern = "acc_", replacement = "")
  
  # check that values in huc12 covars is within range of gage covars
  covar_check <- sw_within_range(gage_covariates,all_huc12_covariates) %>%
    filter(!variable %in% c("comid","huc12","dec_long_va","dec_lat_va","length_km")) %>%
    distinct(index)

  all_huc12_covariates <- all_huc12_covariates[-covar_check$index,]

  return(all_huc12_covariates)
}

sw_gage_all <- function(sites,gage_covariates,fdc_data){
  
  # 1 cfs = 0.028316846592 cms
  cfs2cms <- function(x){x*0.028316846592}
  
  all_gage_data <- fdc_data %>%
    left_join(gage_covariates, by = c("site_no","decade")) %>%
    select(comid,site_no,huc12,decade,dec_long_va,dec_lat_va,everything()) %>%
    arrange(site_no,decade) %>%
    mutate(decade = as.character(decade)) %>%
    filter(decade != 2010) %>%
    na.omit() %>%
    mutate_at(vars(min:L2,median_nonzero), funs(cfs2cms))
}



