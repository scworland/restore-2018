# On 2019-03-11, WHA and ECO discovered that about 60 HUC12 numbers of about
# 9,200 HUC12s in our original releases were no longer existing in a 2017
# USGS Watershed Boundary Dataset (WBD) but were present in the 2016 or
# earlier. The purpose of this script is to read all of our .feather files
# of data and replace the huc12 column with revised values.

library(feather)
library(dataRetrieval)

if(! dir.exists("birdskins/")) dir.create("birdskins/")

HCOV <- read_feather("../../../data/huc12/all_huc12_covariates.feather")
HSOL <- read_feather("../../../data/huc12/all_huc12_solar.feather")

GCOV <- read_feather("../../../data/gage/all_gage_covariates.feather")
GSOL <- read_feather("../../../data/gage/all_gage_solar.feather")
GDAT <- read_feather("../../../data/gage/all_gage_data.feather")
GFLO <- read_feather("../../../data/gage/all_gage_flow_stats.feather")

length(HCOV$comid)         # [1] 55320  (2013-03-12)
length(HSOL$comid)         # [1] 55320  (2013-03-12)
length(unique(HCOV$comid)) # [1]  8988  (2013-03-12)
length(unique(HSOL$comid)) # [1]  8988  (2013-03-12)
length(unique(HCOV$huc12)) # [1]  9220  (2013-03-12)
length(unique(HSOL$huc12)) # [1]  9220  (2013-03-12)

length(GCOV$comid) # [1] 5736  (2013-03-12)
length(GSOL$comid) # [1] 5736  (2013-03-13)
length(GDAT$comid) # [1] 2804  (2013-03-13)
length(GFLO$comid) # [1] 2804  (2013-03-12)

all_comids <-               HCOV$comid
all_comids <- c(all_comids, HSOL$comid)
all_comids <- c(all_comids, GCOV$comid)
all_comids <- c(all_comids, GDAT$comid)
all_comids <- c(all_comids, GFLO$comid)
all_comids <- c(all_comids, GSOL$comid)

all_comids <- unique(all_comids)
length(all_comids) # [1] 9707  (2013-03-12)


# This CSV comes from ECO joining our earlier *huc12*.feather
# to the 2017 WBD.
unzip("201903_all_huc12_covariates.csv.zip", overwrite=TRUE)
WBDh <- read.csv("201903_all_huc12_covariates.csv",
                 header=TRUE, colClasses="character")
unlink("201903_all_huc12_covariates.csv")
head(gsub("\\.0+$", "", WBDh$huc12))
dels <- as.numeric(WBDh$huc12) - as.numeric(WBDh$HUC12_1)
length(dels[abs(dels) > 1])/6
# [1]  198  ECO and WHA spot checked the really big changes. All associated
# it seems with low slopes or drainage basin change of decision.
summary(dels)
quantile(dels, probs=c(0.01,0.99))


length(       HCOV$huc12)    # [1] 55320   (2019-03-13)
length(unique(HCOV$huc12))   # [1]  9220   (2019-03-13)
length(       WBDh$HUC12_1)  # [1] 55320   (2019-03-13)
length(unique(WBDh$HUC12_1)) # [1]  9195   (2019-03-13)

A <- data.frame(huc12=WBDh$HUC12_1, stringsAsFactors=FALSE)
H <- aggregate(A, by=list(A$huc12), length)
J <- H[H$huc12 > 6, ]

# These summary end points (min and maxes) align with the
# parallel code above for the gage statistics
summary(as.numeric(HCOV$huc12) - as.numeric(WBDh$HUC12_1))
summary(as.numeric(HSOL$huc12) - as.numeric(WBDh$HUC12_1))
HCOV$huc12 <- WBDh$HUC12_1
HSOL$huc12 <- WBDh$HUC12_1

write_feather(HCOV, "birdskins/all_huc12_covariates.feather")
write_feather(HSOL, "birdskins/all_huc12_solar.feather")


# These CSVs come from ECO joining our earlier *gage*.feather
# to the 2017 WBD.
unzip("201903_all_gage_covariates.csv.zip")
WBDg <- read.csv("201903_all_gage_covariates.csv",
                 header=TRUE, colClasses="character")
unlink("201903_all_gage_covariates.csv")
length(       GCOV$huc12)    # [1] 5736    (2019-03-13)
length(unique(GCOV$huc12))   # [1]  879    (2019-03-13)
length(       GSOL$huc12)    # [1] 5736    (2019-03-13)
length(unique(GSOL$huc12))   # [1]  879    (2019-03-13)
length(       WBDg$HUC12_1)  # [1] 5736    (2019-03-13)
length(unique(WBDg$HUC12_1)) # [1]  880    (2019-03-13)

dels <- as.numeric(WBDg$huc12) - as.numeric(WBDg$HUC12_1)
length(dels[abs(dels) > 1])/6 # [1]  96 (can divide by 6 like HUC12s)
summary(dels)
quantile(dels, probs=c(0.01,0.99))

# These summary end points (min and maxes) align with the
# parallel code below for the flow statistics
summary(as.numeric(GCOV$huc12) - as.numeric(WBDg$HUC12_1))
summary(as.numeric(GSOL$huc12) - as.numeric(WBDg$HUC12_1))
GCOV$huc12 <- WBDg$HUC12_1
GSOL$huc12 <- WBDg$HUC12_1

unzip("201903_all_gage_flow_stats.csv.zip")
WBDf <- read.csv("201903_all_gage_flow_stats.csv",
                 header=TRUE, colClasses="character")
unlink("201903_all_gage_flow_stats.csv")
length(       GDAT$huc12)    # [1] 2804    (2019-03-13)
length(unique(GDAT$huc12))   # [1]  879    (2019-03-13)
length(       GFLO$huc12)    # [1] 2804    (2019-03-13)
length(unique(GFLO$huc12))   # [1]  879    (2019-03-13)
length(       WBDf$HUC12_1)  # [1] 2804    (2019-03-13)
length(unique(WBDf$HUC12_1)) # [1]  880    (2019-03-13)

dels <- as.numeric(WBDf$huc12) - as.numeric(WBDf$HUC12_1)
length(dels[abs(dels) > 1]) # [1] 65 (can not divide by 6 like HUC12s)
summary(dels)
quantile(dels, probs=c(0.01,0.99))

same_gages    <- unique(GFLO$site_no[GDAT$huc12 == WBDf$HUC12_1])
changed_gages <- unique(GFLO$site_no[GDAT$huc12 != WBDf$HUC12_1])
SG <- readNWISsite(same_gages)
CG <- readNWISsite(changed_gages)
plot(  SG$dec_long_va, SG$dec_lat_va, lwd=0.6, col=8)
points(CG$dec_long_va, CG$dec_lat_va, col=2)

# HUC8 for the Texas site is okay, so no need to tell TXWSC NWIS DBA.
GFLO$huc12[GFLO$site_no == "08075000"]
  CG$huc_cd[ CG$site_no == "08075000"]

  length(same_gages)    # [1] 935
  length(changed_gages) # [1] 21
message("Gage percent change ", round(100*length(changed_gages) /
        (length(same_gages) + length(changed_gages)),digits=2))

# These summary end points (min and maxes) align with the
# parallel code above for the gage statistics
summary(as.numeric(GDAT$huc12) - as.numeric(WBDf$HUC12_1))
summary(as.numeric(GFLO$huc12) - as.numeric(WBDf$HUC12_1))
GDAT$huc12 <- WBDf$HUC12_1
GFLO$huc12 <- WBDf$HUC12_1


write_feather(GCOV, "birdskins/all_gage_covariates.feather")
write_feather(GSOL, "birdskins/all_gage_solar.feather")
write_feather(GDAT, "birdskins/all_gage_data.feather")
write_feather(GFLO, "birdskins/all_gage_flow_stats.feather")

if(dir.exists("./__MACOSX")) unlink("./__MACOSX", recursive=TRUE)
