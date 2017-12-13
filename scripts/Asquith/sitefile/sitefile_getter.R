library(akqdecay); library(sp)

file <- file.choose() # needs definition from ~restore-2018/data/gage
sites <- read.table(file, header=TRUE, colClasses = "character")
sites <- sites$site_no

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

spSites$instruments_cd  <- NULL
spSites$construction_dt <- NULL
spSites$gw_file_cd    <- NULL
spSites$nat_aqfr_cd   <- NULL
spSites$aqfr_cd       <- NULL
spSites$aqfr_type_cd  <- NULL
spSites$well_depth_va <- NULL
spSites$hole_depth_va <- NULL
spSites$depth_src_cd  <- NULL
spSites$project_no    <- NULL

spSites <- spSites[spSites$agency_cd != "USCE",]
# REMOVE "08058900" "USCE" "USCOE E Fk Trinity Rv at McKinney, TX" "ST" 331440 963630 33.2445597 -96.6086019 "M" "F" "NAD27" "NAD83" "48" "48" "085" "US" NA "McKinney East, TX" "  24000" 528.74 "L" 0.1 "NGVD29" "12030106" NA NA "20000907" 164 164 "CST" "Y" NA 164 12783 99 1975 2010 FALSE -56345.0736148357 1130772.79441884 TRUE
# "08058900" "USGS" "E Fk Trinity Rv at McKinney, TX" "ST" 331438 963631 33.24400417 -96.6088797 "M" "F" "NAD27" "NAD83" "48" "48" "085" "US" NA "McKinney East, TX" "  24000" 528.74 "L" 0.1 "NGVD29" "12030106" NA NA NA 164 164 "CST" "Y" NA 164 12783 99 1975 2010 FALSE -56371.1896776191 1130710.93950719 TRUE
# REMOVE "08059400" "USCE" "USCOE Sister Grove Ck nr Blue Ridge, TX" "ST" 331740 962900 33.29455787 -96.4835969 "M" "F" "NAD27" "NAD83" "48" "48" "085" "US" NA "Blue Ridge, TX" "  24000" 526.29 "L" 0.1 "NGVD29" "12030106" NA NA "20000907" 83.1 83.1 "CST" "Y" NA 83.1 10484 70 1975 2016 TRUE -44743.6719397725 1136288.36027692 TRUE
# "08059400" "USGS" "Sister Grove Ck nr Blue Ridge, TX" "ST" 331740 962858 33.29455786 -96.4830413 "M" "F" "NAD27" "NAD83" "48" "48" "085" "US" NA "Blue Ridge, TX" "  24000" 526.29 "L" 0.1 "NGVD29" "12030106" NA NA NA 83.1 83.1 "CST" "Y" NA 83.1 10484 70 1975 2016 TRUE -44692.2667980636 1136288.09772228 TRUE


spSites <- merge(spSites, ExtraSiteFileStuff)

save(spSites, file="spSites.RData")

outfile <- "sitefile.csv"
write.table(spSites, file=outfile, quote=TRUE, row.names=FALSE)

load("spSites.RData")

