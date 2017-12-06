library(akqdecay); library(sp)
load("DV.RData")

ExtraSiteFileStuff <- read.table("ExtraSiteFileStuff.txt", sep="\t", header=TRUE,
                    colClasses="character")
ExtraSiteFileStuff$count_dvs      <- as.numeric(ExtraSiteFileStuff$count_dvs)
ExtraSiteFileStuff$approx_por_pct <- as.numeric(ExtraSiteFileStuff$approx_por_pct)
ExtraSiteFileStuff$begin_cal_year <- as.numeric(ExtraSiteFileStuff$begin_cal_year)
ExtraSiteFileStuff$end_cal_year   <- as.numeric(ExtraSiteFileStuff$end_cal_year)
ExtraSiteFileStuff$active <- as.logical(ExtraSiteFileStuff$active)

ALBEA <- paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 ",
     "+x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
spSites <- sites_to_SpatialPointsDataFrame(sites, proj4string=ALBEA)
spSites$CDA <- spSites$contrib_drain_area_va
spSites$CDA[is.na(spSites$CDA)] <- spSites$drain_area_va[is.na(spSites$CDA)]
spSites$CDA[spSites$CDA == 0] <- spSites$drain_area_va[spSites$CDA == 0]

spSites$gw_file_cd    <- NULL
spSites$nat_aqfr_cd   <- NULL
spSites$aqfr_cd       <- NULL
spSites$aqfr_type_cd  <- NULL
spSites$well_depth_va <- NULL
spSites$hole_depth_va <- NULL
spSites$depth_src_cd  <- NULL
spSites$project_no    <- NULL

spSites <- merge(spSites, ExtraSiteFileStuff)

save(spSites, file="spSites.RData")

write.table(spSites, file="sitefile.csv", quote=TRUE, row.names=FALSE)

load("spSites.RData")

