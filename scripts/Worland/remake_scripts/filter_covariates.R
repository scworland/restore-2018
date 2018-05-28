
df <- remake::fetch('gage_all')

sw_gage_filter <- function(df) {
  
  # used to convert acre-ft to km^2-m
  val <- 0.3048*(1/247.104393047)
  
  # list of sites impacted by edwards aquifer
  edwards_sites <- c("08156800","08155300","08155400","08181400","08185000",
                     "08190500","08197500","08198500","08200700","08202700")
  
  # 1 cfs = 0.028316846592 cms
  cfs2cms <- function(x){x*0.028316846592}
  
  # edwards_comids <- c("5781337","5781711","5781703","10835030","7851041",
  #                     "7876116","10646917","10646003","10655659","10655687")

  
  df_clean <- df %>%
    mutate(flood_storage = abs((acc_nid_storage*val - acc_norm_storage*val)/acc_basin_area),
           aquifers = ifelse(aquifers == "cat_aq_nodata","nodata",aquifers),
           soller = ifelse(soller == "cat_soller_nodata","nodata",soller),
           ed_rch_zone = ifelse(site_no %in% edwards_sites, "1", "0")) %>%
    select(-acc_elev_max,-acc_elev_min,-runoff_mean,-runoff_sd,-acc_hdens,
           -acc_nid_storage,-acc_norm_storage) %>%
    mutate_at(vars(min:L2), funs(cfs2cms))
    
  
  
}

sw_huc12_filter <- function(df) {
  
}

gages <- remake::fetch('gage_covariates')
hucs <- remake::fetch('huc12_covariates')

hold <- sw_within_range(gages,hucs) %>%
  filter(!variable %in% c("comid","huc12","lon","lat"))






