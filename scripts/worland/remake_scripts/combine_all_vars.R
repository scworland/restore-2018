
sw_combine_gage <- function(sites,immutable_vars,dams,housing_density,climate,lulc){
  
  filter <- dplyr::filter
  select <- dplyr::select
  
  ## variables to be combined
  immutable_gage_df <- immutable_vars$immutable_gage_df 
  house_gage_df <- housing_density$house_gage_df
  dam_gage_df <- dams$dam_gage_df 
  climate_gage_df <- climate$climate_gage_df
  lulc_gage_df <- lulc$lulc_gage_df
  
  all_gage_covariates <- list(house_gage_df,dam_gage_df,climate_gage_df,lulc_gage_df) %>%
    reduce(left_join, by = c("comid","site_no","decade")) %>%
    arrange(comid,decade) %>%
    left_join(immutable_gage_df, by = c("comid","site_no")) %>%
    left_join(sites,by=c("comid","site_no")) %>%
    rename_all(tolower) %>%
    select(comid,site_no,huc12,lon,lat,acc_hdens:woody_wetland,
           acc_bfi:acc_twi,acc_basin_area:acc_rdx,everything()) %>%
    filter(decade %in% c(1950,1960,1970,1980,1990,2000))
  
  return(all_gage_covariates)
}

sw_combine_huc12 <- function(huc12s,immutable_vars,dams,housing_density,climate,lulc){
  
  ## variables to be combined
  immutable_huc12_df <- immutable_vars$immutable_huc12_df 
  house_huc12_df <- housing_density$house_huc12_df
  dam_huc12_df <- dams$dam_huc12_df 
  climate_huc12_df <- climate$climate_huc12_df
  lulc_huc12_df <- lulc$lulc_huc12_df
  
  all_huc12_covariates <- list(house_huc12_df,dam_huc12_df,climate_huc12_df,lulc_huc12_df) %>%
    reduce(left_join, by = c("comid","huc12","decade")) %>%
    arrange(comid,decade) %>%
    left_join(immutable_huc12_df, by = c("comid","huc12")) %>%
    left_join(huc12s, by = c("comid","huc12")) %>%
    rename_all(tolower) %>%
    select(comid,huc12,lon,lat,acc_hdens:woody_wetland,
           acc_bfi:acc_twi,acc_basin_area:acc_rdx,everything()) %>%
    filter(decade %in% c(1950,1960,1970,1980,1990,2000))

  return(all_huc12_covariates)
}

sw_gage_all <- function(sites,gage_covariates,fdc_data){
  all_gage_data <- fdc_data %>%
    left_join(gage_covariates, by = c("site_no","decade")) %>%
    select(comid,site_no,huc12,decade,lon,lat,everything()) %>%
    arrange(site_no,decade) %>%
    mutate(decade = as.character(decade)) %>%
    filter(decade != 2010) %>%
    na.omit()
}



