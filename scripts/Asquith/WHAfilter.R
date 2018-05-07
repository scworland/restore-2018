gage_covariates_WHAfilter <- function(df) {

  if(any(is.na(df$comid))) { # SCW has already done this?
     stop("COMID failure")
  }
  #length(df$comid)                           # [1] 2804
  #length(df$comid[is.na(df$comid)])          # [1] 0
  df <- df[! is.na(df$comid),] # SCW has already done this?
  #length(df$comid[is.na(df$acc_basin_area)]) # [1] 0
  df <- df[! is.na(df$acc_basin_area), ] # SCW has already done this?
  #length(df$comid[is.na(df$area_sqkm)])      # [1] 0
  df <- df[! is.na(df$area_sqkm), ] # SCW has already done this? What is the meaning with this?

  # *** ABOUT SLOPE ***
  # df[df$acc_basin_slope <= 0,] # We appear to have none but we have some
  # zeros for prediction COMIDs. These we will reset these to a lower bounds of 0.0005.
  #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 0.020   1.940   3.480   4.836   6.130  34.080
  # We made need to revisit a lower limit for the prediction COMIDs.

  # *** ABOUT WHA's "FLOOD STORAGE" ***
  # storages are in acre-ft
  # 1 km2 = 247.104393047 acres
  df$acc_nid_storage <- df$acc_nid_storage*0.3048*(1/247.104393047) # to km^2-meter storage
  df$acc_norm_storage <- df$acc_norm_storage*0.3048*(1/247.104393047) # to km^2-meter storage
  df$flood_storage <- (df$acc_nid_storage - df$acc_norm_storage) / df$acc_basin_area # meters of storage

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

  #df$log_flood_storage <- log10(df$flood_storage+.001/0.3048) # Asquith testing only (100th of a foot)
  #plot(qnorm(lmomco::pp(df$log_flood_storage)), sort(df$log_flood_storage)) # Asquith testing only
  #mtext(paste0("Maximum log10offets of flood_storage=",max(df$log_flood_storage))) # Asquith testing only
  # The conclusion is that a maximum log10 flood storage is 0 (1 m, or about 3.3 ft). We will
  # be truncated prediction COMIDs according to this upper limit.
  df <- df[log10(df$flood_storage) < 0,]

  # *** ABOUT PERMANENT ICE AND SNOW ***
  # We have no permanent ice and snow, just remove the column entirely.
  # perennial_ice_snow
  df$perennial_ice_snow <- NULL
  df$aquifers[df$aquifers == "cat_aq_nodata"] <- "nodata"
  df$soller[df$soller == "cat_soller_nodata"] <- "nodata"

  # Checks on existance of the factors we need. Here we list them all.
  message("GAGE_COVARIATES: bedperm length of levels: ", length(levels(as.factor(df$bedperm))))
  message("GAGE_COVARIATES: aquifers length of levels: ",length(levels(as.factor(df$aquifers))))
  message("GAGE_COVARIATES: soller length of levels: ",  length(levels(as.factor(df$soller))))
  message("GAGE_COVARIATES: hlr length of levels: ",     length(levels(as.factor(df$hlr))))
  message("GAGE_COVARIATES: ecol3 length of levels: ",   length(levels(as.factor(df$ecol3))))
  message("GAGE_COVARIATES: physio length of levels: ",  length(levels(as.factor(df$physio))))
  message("GAGE_COVARIATES: statsgo length of levels: ", length(levels(as.factor(df$statsgo))))

  message("GAGE_COVARIATES: bedperm length of levels: ", paste(levels(as.factor(df$bedperm)), collapse=" "))
  message("GAGE_COVARIATES: aquifers length of levels: ",paste(levels(as.factor(df$aquifers)), collapse=" "))
  message("GAGE_COVARIATES: soller length of levels: ",  paste(levels(as.factor(df$soller)), collapse=" "))
  message("GAGE_COVARIATES: hlr length of levels: ",     paste(levels(as.factor(df$hlr)), collapse=" "))
  message("GAGE_COVARIATES: ecol3 length of levels: ",   paste(levels(as.factor(df$ecol3)), collapse=" "))
  message("GAGE_COVARIATES: physio length of levels: ",  paste(levels(as.factor(df$physio)), collapse=" "))
  message("GAGE_COVARIATES: statsgo length of levels: ", paste(levels(as.factor(df$statsgo)), collapse=" "))

  # *** ABOUT TEXAS Edwards Aquifer Hydrology ***  (We never got clarity concerning Florida.)
  # Edwards Recharge Zone in Texas has large affect to the low end of the streamflow
  # regime. WHA after consulting a shapefile, the high no flow percentages, and
  # basic situational awareness of the region, has identified these as needed a flag.
  # We have a shapefile for this. These are not all precisely in the zone as some
  # of the gages are a little downstream as part of regional water balance needs.
  df$edwards_recharge_zone_impacted <- rep(0, length(df$site_no))
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
  df$edwards_recharge_zone_impacted <- as.character(df$edwards_recharge_zone_impacted)

  # Flow-duration curve and statistics to cms (1 cfs = 0.028316846592 cms)
  # WHA is well aware this is long winded but doesn't want to use column index for safety
  # and does not want to convert the data.frame into an environment and list just in case.
  P <- 0.028316846592
  df$min    <- P*df$min
  df$f0.02  <- P*df$f0.02
  df$f0.05  <- P*df$f0.05
  df$f0.1   <- P*df$f0.1
  df$f0.2   <- P*df$f0.2
  df$f0.5   <- P*df$f0.5
  df$f01    <- P*df$f01
  df$f02    <- P*df$f02
  df$f05    <- P*df$f05
  df$f10    <- P*df$f10
  df$f20    <- P*df$f20
  df$f25    <- P*df$f25
  df$f30    <- P*df$f30
  df$f40    <- P*df$f40
  df$f50    <- P*df$f50
  df$f60    <- P*df$f60
  df$f70    <- P*df$f70
  df$f75    <- P*df$f75
  df$f80    <- P*df$f80
  df$f90    <- P*df$f90
  df$f95    <- P*df$f95
  df$f98    <- P*df$f98
  df$f99    <- P*df$f99
  df$f99.5  <- P*df$f99.5
  df$f99.8  <- P*df$f99.8
  df$f99.9  <- P*df$f99.9
  df$f99.95 <- P*df$f99.95
  df$f99.98 <- P*df$f99.98
  df$max    <- P*df$max
  df$L1     <- P*df$L1
  df$L2     <- P*df$L2
  df$median_nonzero <- df$median_nonzero

  return(df)
}

huc12_covariates_WHAfilter <- function(df, purge=FALSE) {
  #length(df$comid[df$decade == 1950])                           # [1] 10302
  #length(df$comid[is.na(df$comid) & df$decade == 1950])          # [1] 258
  df <- df[! is.na(df$comid),] # SCW has already done this?
  #length(df$comid[is.na(df$acc_basin_area) & df$decade == 1950]) # [1] 61
  df <- df[! is.na(df$acc_basin_area), ] # SCW has already done this?
  #length(df$comid[is.na(df$area_sqkm) & df$decade == 1950])      # [1] 36
  df <- df[! is.na(df$area_sqkm), ] # SCW has already done this? What is the meaning with this?

  # *** ABOUT WHA's "FLOOD STORAGE" ***
  # storages are in acre-ft
  # 1 km2 = 247.104393047 acres
  df$acc_nid_storage <- df$acc_nid_storage*0.3048*(1/247.104393047) # to km^2-meter storage
  df$acc_norm_storage <- df$acc_norm_storage*0.3048*(1/247.104393047) # to km^2-meter storage
  df$flood_storage <- (df$acc_nid_storage - df$acc_norm_storage) / df$acc_basin_area # meters of storage
  df$flood_storage[is.na(df$flood_storage)] <- 0 # resetting the missing
  #df[df$flood_storage < 0,]; # plot(df); plot(df[df$flood_storage < 0,], add=TRUE, col=2)
  #summary(df$flood_storage[df$flood_storage < 0])
  #  Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
  #-3.53830 -0.08787 -0.05581 -0.50040 -0.04884 -0.02663
  #length(df$flood_storage[df$flood_storage < 0]) # [1] 10
  df$flood_storage <- abs(df$flood_storage) # there appear to be spurious negatives, switched values in original data?
  #df$log_flood_storage <- log10(df$flood_storage+.001/0.3048) # Asquith testing only (100th of a foot)
  #plot(qnorm(lmomco::pp(df$log_flood_storage)), sort(df$log_flood_storage), type="l") # Asquith testing only
  #lines(qnorm(lmomco::pp(DD$log_flood_storage)), sort(DD$log_flood_storage), col=2) # Asquith testing only
  #length(df$comid[df$log_flood_storage >= 0]) # This is 1 meter or about 3.3 feet of watershed depth storage
  # [1] 35
  #summary(df$flood_storage[df$flood_storage >= 0]) # Asquith testing only
  #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 0.5228  0.5546  0.7193  0.7827  0.8425  1.5351
  df <- df[log10(df$flood_storage) < 0,] # COMID parameter space restriction

  # *** ABOUT A MINIMUM SLOPE ***
  # WHA worked on a TxDOT project many years ago to identify a minimum slope for statistical
  # operations. What happens as slope ---> zero is that gravitational forces go away and
  # pressure differentials are the cause of the flow. We could mitigate for this decision
  # here by making a truncation. The TxDOT work used a minimum of 0.0005 dimensionless, so
  # we scale here to percent. However, I am explicitly using 0.0002 as this is the minimum
  # of our gages and pretty close to 0.0005 anyway.
  df$acc_basin_slope[df$acc_basin_slope == 0] <- 0.0002*100 # COMID parameter space restriction
  # One open problem is that I am showing 1,788 COMIDs with a missing slope. Crap what should we do?
  # RECOMMENDATION: We need a lower slope trunction on the COMIDs in my mind but not on the gages
  # to build the model.

  # *** ABOUT PERMANENT ICE AND SNOW ***
  # We have no permanent ice and snow, just remove the column entirely.
  # perennial_ice_snow
  df$perennial_ice_snow <- NULL
  df$aquifers[df$aquifers == "cat_aq_nodata"] <- "nodata"
  df$soller[df$soller == "cat_soller_nodata"] <- "nodata"

  # ***** WHA/RRK mental notes and checking the ungaged categories
  #CAT_SOLLER_101 ---> CAT_SOLLER_11   *** No we could delete, only one? ***
  #CAT_SOLLER_11 Estimated percent of catchment that contains the surficial materials: Alluvial sediments, thin, less than 100 feet thick; Holocene to Pliocene geologic age.
  #to
  #CAT_SOLLER_101 Estimated percent of catchment that contains the surficial materials:Biological sediments less than 100 feet thick. Holocene to middle Pleistocene geologic age.
  # COV[! is.na(COV$hlr) & COV$soller == "cat_soller_101" & COV$decade == 1950,2] #  A tibble: 1 x 1

  #CAT_HLR_17 ---> CAT_HLR_* 15/18
  #CAT_HLR_ * 15/18
  #to
  #CAT_HLR_17
  #COV[! is.na(COV$hlr) & COV$hlr == "cat_hlr_17" & COV$decade == 1950,2] # A tibble: 8 x 1

  # Checks on existance of the factors we need. Here we list them all.
  message("HUC12_COVARIATES: bedperm length of levels: ", length(levels(as.factor(df$bedperm))))
  message("HUC12_COVARIATES: aquifers length of levels: ",length(levels(as.factor(df$aquifers))))
  message("HUC12_COVARIATES: soller length of levels: ",  length(levels(as.factor(df$soller))))
  message("HUC12_COVARIATES: hlr length of levels: ",     length(levels(as.factor(df$hlr))))
  message("HUC12_COVARIATES: ecol3 length of levels: ",   length(levels(as.factor(df$ecol3))))
  message("HUC12_COVARIATES: physio length of levels: ",  length(levels(as.factor(df$physio))))
  message("HUC12_COVARIATES: statsgo length of levels: ", length(levels(as.factor(df$statsgo))))

  message("HUC12_COVARIATES: bedperm length of levels: ", paste(levels(as.factor(df$bedperm)), collapse=" "))
  message("HUC12_COVARIATES: aquifers length of levels: ",paste(levels(as.factor(df$aquifers)), collapse=" "))
  message("HUC12_COVARIATES: soller length of levels: ",  paste(levels(as.factor(df$soller)), collapse=" "))
  message("HUC12_COVARIATES: hlr length of levels: ",     paste(levels(as.factor(df$hlr)), collapse=" "))
  message("HUC12_COVARIATES: ecol3 length of levels: ",   paste(levels(as.factor(df$ecol3)), collapse=" "))
  message("HUC12_COVARIATES: physio length of levels: ",  paste(levels(as.factor(df$physio)), collapse=" "))
  message("HUC12_COVARIATES: statsgo length of levels: ", paste(levels(as.factor(df$statsgo)), collapse=" "))
  if(! purge) return(df)

  # Never gaged, we need to remove these categories from the corresponding variable.
  # bedperm  *** NO ACTION TO BE DONE ***
  # aquifers cat_aq307 cat_aq405 cat_aq407
  # sroller cat_soller_101 cat_soller_321 cat_soller_322 cat_soller_811 cat_soller_820 cat_soller_822
  # hlr cat_hlr_17 cat_hlr_19
  # ecol3 ecol3_37 ecol3_39  ecol3_71 ecol3_72 nodata
  # physio cat_physio_11 cat_physio_14 cat_physio_5
  # statsgo cat_hgad

  length(df$comid[df$aquifers == "cat_aq307" & df$decade == 1950]) # [1] 1
  length(df$comid[df$aquifers == "cat_aq405" & df$decade == 1950]) # [1] 19
  length(df$comid[df$aquifers == "cat_aq407" & df$decade == 1950]) # [1] 6
  df <- df[df$aquifers != "cat_aq307", ]
  df <- df[df$aquifers != "cat_aq405", ]
  df <- df[df$aquifers != "cat_aq407", ]

  length(df$comid[df$soller == "cat_soller_101" & df$decade == 1950]) # [1] 1
  length(df$comid[df$soller == "cat_soller_321" & df$decade == 1950]) # [1] 4
  length(df$comid[df$soller == "cat_soller_322" & df$decade == 1950]) # [1] 1
  length(df$comid[df$soller == "cat_soller_811" & df$decade == 1950]) # [1] 2
  length(df$comid[df$soller == "cat_soller_820" & df$decade == 1950]) # [1] 69
  length(df$comid[df$soller == "cat_soller_822" & df$decade == 1950]) # [1] 1
  df <- df[df$soller != "cat_soller_101", ]
  df <- df[df$soller != "cat_soller_321", ]
  df <- df[df$soller != "cat_soller_322", ]
  df <- df[df$soller != "cat_soller_811", ]
  df <- df[df$soller != "cat_soller_820", ]
  df <- df[df$soller != "cat_soller_822", ]

  length(df$comid[df$hlr == "cat_hlr_17" & df$decade == 1950]) # [1] 8
  length(df$comid[df$hlr == "cat_hlr_19" & df$decade == 1950]) # [1] 4
  df <- df[df$hlr != "cat_hlr_17", ]
  df <- df[df$hlr != "cat_hlr_19", ]

  length(df$comid[df$ecol3 == "ecol3_37" & df$decade == 1950]) # [1] 51
  length(df$comid[df$ecol3 == "ecol3_39" & df$decade == 1950]) # [1] 5
  length(df$comid[df$ecol3 == "ecol3_71" & df$decade == 1950]) # [1] 8
  length(df$comid[df$ecol3 == "ecol3_72" & df$decade == 1950]) # [1] 7
  length(df$comid[df$ecol3 == "nodata" & df$decade == 1950]) # [1] 13
  df <- df[df$ecol3 != "ecol3_37", ]
  df <- df[df$ecol3 != "ecol3_39", ]
  df <- df[df$ecol3 != "ecol3_71", ]
  df <- df[df$ecol3 != "ecol3_72", ]
  df <- df[df$ecol3 != "nodata", ]

  length(df$comid[df$physio == "cat_physio_11" & df$decade == 1950]) # [1] 5
  length(df$comid[df$physio == "cat_physio_14" & df$decade == 1950]) # [1] 3
  length(df$comid[df$physio == "cat_physio_5" & df$decade == 1950]) # [1] 38
  df <- df[df$physio != "cat_physio_11", ]
  df <- df[df$physio != "cat_physio_14", ]
  df <- df[df$physio != "cat_physio_5", ]

  length(df$comid[df$statsgo == "cat_hgad" & df$decade == 1950]) # [1] 8
  df <- df[df$statsgo != "cat_hgad", ]

  message(" ******************************************* ")
  # Checks on existance of the factors we need. Here we list them all.
  message("HUC12_COVARIATES: bedperm length of levels: ", length(levels(as.factor(df$bedperm))))
  message("HUC12_COVARIATES: aquifers length of levels: ",length(levels(as.factor(df$aquifers))))
  message("HUC12_COVARIATES: soller length of levels: ",  length(levels(as.factor(df$soller))))
  message("HUC12_COVARIATES: hlr length of levels: ",     length(levels(as.factor(df$hlr))))
  message("HUC12_COVARIATES: ecol3 length of levels: ",   length(levels(as.factor(df$ecol3))))
  message("HUC12_COVARIATES: physio length of levels: ",  length(levels(as.factor(df$physio))))
  message("HUC12_COVARIATES: statsgo length of levels: ", length(levels(as.factor(df$statsgo))))

  message("HUC12_COVARIATES: bedperm length of levels: ", paste(levels(as.factor(df$bedperm)), collapse=" "))
  message("HUC12_COVARIATES: aquifers length of levels: ",paste(levels(as.factor(df$aquifers)), collapse=" "))
  message("HUC12_COVARIATES: soller length of levels: ",  paste(levels(as.factor(df$soller)), collapse=" "))
  message("HUC12_COVARIATES: hlr length of levels: ",     paste(levels(as.factor(df$hlr)), collapse=" "))
  message("HUC12_COVARIATES: ecol3 length of levels: ",   paste(levels(as.factor(df$ecol3)), collapse=" "))
  message("HUC12_COVARIATES: physio length of levels: ",  paste(levels(as.factor(df$physio)), collapse=" "))
  message("HUC12_COVARIATES: statsgo length of levels: ", paste(levels(as.factor(df$statsgo)), collapse=" "))

  return(df)
}


