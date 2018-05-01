gage_covariates_WHAfilter <- function(df) {
  # *** ABOUT SLOPE ***
  # df[df$acc_basin_slope <= 0,] # We appear to have none but we have some
  # zeros for prediction COMIDs. These we will reset these to a lower bounds of 0.0005.
  #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 0.020   1.940   3.480   4.836   6.130  34.080
  # We made need to revisit a lower limit for the prediction COMIDs.

  # *** ABOUT WHA's "FLOOD STORAGE" ***
  # storages are in acre-ft
  # 1 km2 = 247.104393047 acres
  df$flood_storage <- (df$acc_nid_storage - df$acc_norm_storage) /
                              (df$acc_basin_area*247.104393047)

  # Two streamgages both had one decade (2000) for which the flood storage is a negative
  # value, and these streamgages are in Florida. This implies that an inconsistency in
  # source USACE NID data or processing into the products of \citep{NHDplus} including
  # streamgage location on the virtual network.

  # For streamgage 02295420, the computed flood storage from the raw data was
  # $-3.344$~feet but the value for 1990 was $0.400$~feet and for 1980 was $0.275$~feet.
  # For this study, the value for 2000 was assigned that of 1990; held constant, in other
  # words. For streamgage 02296750, the computed flood storage from the raw data was
  # $-0.174$~feet and for 1980 was $0.142$~feet and for 1990 was $0.158$~feet. For this
  # study, the value for 2000 was assigned as the absolute value for an assumption that
  # data transcription error in the NID is involved.
  # df[DD$site_no == "02295420",] # Asquith testing only
  # df[DD$site_no == "02296750",] # Asquith testing only
  df$flood_storage[                 df$site_no == "02295420" & df$decade == "2000"] <-
                   df$flood_storage[df$site_no == "02295420" & df$decade == "1990"]
  df$flood_storage[                 df$site_no == "02296750" & df$decade == "2000"] <-
               abs(df$flood_storage[df$site_no == "02296750" & df$decade == "2000"])

  #df$flood_storage <- log10(df$flood_storage+.001) # Asquith testing only
  #plot(qnorm(lmomco::pp(df$flood_storage)), sort(df$flood_storage)) # Asquith testing only
  #text(paste0("Maximum log10offets of flood_storage=",max(df$flood_storage))) # Asquith testing only
  # The conclusion is that a maximum log10 flood storage is 0.52 (3.3 ft). We will
  # be truncated prediction COMIDs according to this upper limit.

  # *** ABOUT PERMANENT ICE AND SNOW ***
  # We have no permanent ice and snow, just remove the column entirely.
  # perennial_ice_snow
  df$perennial_ice_snow <- NULL

  # *** ABOUT TEXAS Edwards Aquifer Hydrology ***  (We never got clarity about Florida.)
  # Edwards Recharge Zone in Texas has large affect to the low end of the streamflow
  # regime. WHA after consulting a shapefile, the high no flow percentages, and
  # basic situational awareness of the region, has identified these as needed a flag.
  # We have a shapefile for this. These are not all precisely in the zone as some
  # of the gages are a little downstream as part of regional water balance needs.
  df$edwards_recharge_zone_impacted <- 0
  df$edwards_recharge_zone_impacted[df$site_no == "08156800"] <- 1
  df$edwards_recharge_zone_impacted[df$site_no == "08155300"] <- 1
  df$edwards_recharge_zone_impacted[df$site_no == "08155400"] <- 1
  df$edwards_recharge_zone_impacted[df$site_no == "08181400"] <- 1
  df$edwards_recharge_zone_impacted[df$site_no == "08184000"] <- 1
  df$edwards_recharge_zone_impacted[df$site_no == "08185000"] <- 1
  df$edwards_recharge_zone_impacted[df$site_no == "08190500"] <- 1
  df$edwards_recharge_zone_impacted[df$site_no == "08197500"] <- 1
  df$edwards_recharge_zone_impacted[df$site_no == "08198500"] <- 1
  df$edwards_recharge_zone_impacted[df$site_no == "08200700"] <- 1
  df$edwards_recharge_zone_impacted[df$site_no == "08202700"] <- 1
  df$edwards_recharge_zone_impacted <- as.factor(df$edwards_recharge_zone_impacted)

  return(df)
}

huc12_covariates_WHAfilter <- function(df) {
  #length(df$comid)                           # [1] 59070
  #length(df$comid[is.na(df$comid)])          # [1] 1458
  df <- df[! is.na(df$comid),] # SCW has already done this?
  #length(df$comid[is.na(df$acc_basin_area)]) # [1] 330
  df <- df[! is.na(df$acc_basin_area), ] # SCW has already done this?
  #length(df$comid[is.na(df$area_sqkm)])      # [1] 210
  df <- df[! is.na(df$area_sqkm), ] # SCW has already done this? What is the meaning with this?

  # *** ABOUT WHA's "FLOOD STORAGE" ***
  # storages are in acre-ft
  # 1 km2 = 247.104393047 acres
  df$flood_storage <- df$acc_nid_storage - df$acc_norm_storage
  df$flood_storage <- df$flood_storage/(df$acc_basin_area*247.104393047)
  df$flood_storage[is.na(df$flood_storage)] <- 0 # resetting the missing
  #df[df$flood_storage < 0,]; # plot(df); plot(df[df$flood_storage < 0,], add=TRUE, col=2)
  df$flood_storage <- abs(df) # there appear to be spurious negatives, switched values in original data?
  #df$flood_storage <- log10(df$flood_storage+.01) # Asquith testing only
  #plot(qnorm(pp(df$flood_storage)), sort(df$flood_storage), type="l") # Asquith testing only
  #lines(qnorm(pp(DD$flood_storage)), sort(DD$flood_storage), col=2) # Asquith testing only
  #length(df$comid[df$flood_storage >= 0.52]) # This is about 3.3 feet of watershed depth storage
  # [1] 35
  #summary(df$flood_storage[df$flood_storage >= 0.52]) # Asquith testing only
  #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 0.5228  0.5546  0.7193  0.7827  0.8425  1.5351
  df <- df[df$flood_storage < 0.52,] # COMID parameter space restriction

  # *** ABOUT A MINIMUM SLOPE ***
  # WHA worked on a TxDOT project many years ago to identify a minimum slope for statistical
  # operations. What happens as slope ---> zero is that gravitational forces go away and
  # pressure differentials are the cause of the flow. We could mitigate for this decision
  # here by making a truncation. The TxDOT work used a minimum of 0.0005 dimensionless, so
  # we scale here to percent. However, I am explicitly using 0.0002 as this is the minimum
  # of our gages and pretty close to 0.0005.
  df$acc_basin_slope[df$acc_basin_slope == 0] <- 0.0002*100 # COMID parameter space restriction
  # One open problem is that I am showing 1,788 COMIDs with a missing slope. Crap what should we do?
  # RECOMMENDATION: We need a lower slope trunction on the COMIDs in my mind but not on the gages
  # to build the model.

  # WHA's alpha prediction testing fails at the predict.gam() level because some of these
  # factors in the prediction COMIDs have the category but the observational dataset
  # does not have anything. The model will know nothing about this categories in otherwords.
  #message("REMOVING nodata (Bed Permeability)")
  #length(df$comid[df$bedperm == "nodata"]) # [1] 198
  df <- df[df$bedperm != "nodata",]

  # WHA's alpha prediction testing fails at the predict.gam() level because some of these
  # factors in the prediction COMIDs have the category but the observational dataset
  # does not have anything. The model will know nothing about this categories in otherwords.
  #message("REMOVING ecol3_37, ecol3_72, nodata (Ecoregion)")
  #length(df$comid[df$ecol3 == "ecol3_37"]) # [1] 198
  #length(df$comid[df$ecol3 == "ecol3_72"]) # [1] 12
  #length(df$comid[df$ecol3 == "nodata"])   # [1] 30
  df <- df[df$ecol3 != "ecol3_37",]
  df <- df[df$ecol3 != "ecol3_72",]
  df <- df[df$ecol3 != "nodata",  ]

  return(df)
}


