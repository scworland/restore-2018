library(akqdecay); library(sp)

ALBEA <- paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 ",
     "+x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
spModSites <- sites_to_SpatialPointsDataFrame(ModelSites, proj4string=ALBEA)
spModSites$CDA <- spModSites$contrib_drain_area_va
spModSites$CDA[is.na(spModSites$CDA)] <- spModSites$drain_area_va[is.na(spModSites$CDA)]
spModSites$CDA[spModSites$CDA == 0] <- spModSites$drain_area_va[spModSites$CDA == 0]

spModSites$instruments_cd  <- NULL
spModSites$construction_dt <- NULL
spModSites$gw_file_cd    <- NULL
spModSites$nat_aqfr_cd   <- NULL
spModSites$aqfr_cd       <- NULL
spModSites$aqfr_type_cd  <- NULL
spModSites$well_depth_va <- NULL
spModSites$hole_depth_va <- NULL
spModSites$depth_src_cd  <- NULL
spModSites$project_no    <- NULL

