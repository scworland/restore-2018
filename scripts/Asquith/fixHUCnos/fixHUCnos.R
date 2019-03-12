# On 2019-03-11, WHA and ECO discovered that about 60 HUC12 numbers of about
# 9,200 HUC12s in our original releases were no longer existing in a 2017
# USGS Watershed Boundary Dataset (WBD) but were present in the 2016 or
# earlier. The purpose of this script is to read all of our .feather files
# of data and replace the huc12 column with revised values.

library(feather)

HCOV <- read_feather("../../../data/huc12/all_huc12_covariates.feather")
HSOL <- read_feather("../../../data/huc12/all_huc12_solar.feather")

GCOV <- read_feather("../../../data/gage/all_gage_covariates.feather")
GDAT <- read_feather("../../../data/gage/all_gage_data.feather")
GFLO <- read_feather("../../../data/gage/all_gage_flow_stats.feather")
GSOL <- read_feather("../../../data/gage/all_gage_solar.feather")

length(HCOV$comid) # [1] 55320  (2013-03-12)
length(HSOL$comid) # [1] 55320  (2013-03-12)

length(GCOV$comid) # [1] 5736  (2013-03-12)
length(GDAT$comid) # [1] 5736  (2013-03-12)
length(GFLO$comid) # [1] 2804  (2013-03-12)
length(GSOL$comid) # [1] 2804  (2013-03-12)

all_comids <-               HCOV$comid
all_comids <- c(all_comids, HSOL$comid)
all_comids <- c(all_comids, GCOV$comid)
all_comids <- c(all_comids, GDAT$comid)
all_comids <- c(all_comids, GFLO$comid)
all_comids <- c(all_comids, GSOL$comid)

all_comids <- unique(all_comids)
length(all_comids) # [1] 9707  (2013-03-12)


