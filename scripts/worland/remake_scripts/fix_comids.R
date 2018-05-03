
sw_fix_comids <- function(site_comids){
  
  st_francis <- c("07035500","07035800","07037500","07037700","07039500")
  
  clean_sites <- site_comids %>%
    mutate(comid = ifelse(site_no=="02301745","16918094",comid),
           comid = ifelse(site_no=="07032200","14199477",comid),
           comid = ifelse(site_no=="07045000","766886",comid),
           comid = ifelse(site_no=="07281000","18050376",comid),
           comid = ifelse(site_no=="07335400","8348979",comid),
           comid = ifelse(site_no=="02326900","10323850",comid),
           comid = ifelse(site_no=="02310800","16954130",comid),
           comid = ifelse(site_no=="07044000","766886",comid),
           comid = ifelse(site_no=="07341500","19947688",comid),
           comid = ifelse(site_no=="02359315","2178743",comid),
           comid = ifelse(site_no=="02301740","16919448",comid)) %>%
    filter(!site_no %in% st_francis) %>% # remove st. francis gages
    na.omit()
  
  return(clean_sites)
}


# Additional notes 5/1/2018
# site number "07341500" is associated with comid "19947688" because the "updated_locs"
# file erroneously dropped the site because someone thought it was "stage only". So
# we added it back in here.
#
# Mike Wieczorek sent email saying site_no "02359315" has comid "2178743" and
# site_no "02301740" has comid "16919448" but drainage area is slightly different
# than NWIS (6 mi^2 vs 2 mi^2)

# These are notes concerning COMID realignment for streamgages having a large descrepancy
# between NWIS reported area and the acc_basin_area from the NHDplus
# (Wieczorek and others, 2017) via the all_gage_data.feather file (from Scott Worland)
# of flow-duration curves and covariates. A check against the drainage area (minimum of
# total and contributing) to the acc_basin_area was made in
# scripts/Asquith/fdc_model/fdcest.R.
# 
# We are recommending revision to COMIDs as per discussion below. We have 11 sites of which
# 3 are not to be changed. We must tentatively conclude that NWIS areas are at "fault" and
# not to be investigated further. We think that eight can be fixed, and these often have
# area differences measured in orders of magnitude.
# 
# A large descrepancy is measured in log10-cycle differences.
# 
# -WHA, -RRK (April 24, 2018)
# 
# REFERENCES
# \bibitem[U.S. Geological Survey(2014)]{StreamsOnly}
# U.S. Geological Survey, 2014. USGS Small-scale dataset---1:1,000,000-scale streams of the
# United States 201403 shapefile: U.S. Geological Survey,
# \url{https://www.sciencebase.gov/catalog/item/581d052de4b08da350d524eb}
# 
# \bibitem[U.S. Geological Survey(2018)]{WBD}
# U.S. Geological Survey, 2018. Watershed boundary dataset, accessed March 30, 2018, at
# \url{https://nhd.usgs.gov/wbd.html}.
# 
# \bibitem[Wieczorek et~al.(2017)]{NHDplus}
# Wieczorek, M.E., Jackson, S.E., and Schwarz, G.E., 2017. Select attributes for NHDPlus
# version 2.1 reach catchments and modified network routed upstream watersheds for the
# conterminous United States, (dated September 30, 2017 [draft]), U.S. Geological Survey
# data release, accessed on March 12, 2018 at \url{https://doi.org/10.5066/F7765D7V}.
# 
# 
# 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   +   NEW SECTION. Outliers of |>1/2| log10 and most are several orders of magnitude off.
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   ===***===***===***===***===***===***===***===***===***===***===***===***===***===***
#   coordinates  site_no    comid        huc12 decade       lon      lat    n
# (1341699, 633499.5) 02301745 16918140 031002060302   2000 -82.37704 27.90225 3653
# 02301745 DELANEY CREEK POPOFF CANAL NEAR TAMPA FL
# point does not land on line work
# most probable COMID based on shifting point due north is 16918094 (Delaney Creek)
# most DS COMID for HUC12: 16918108
# HUC12 agrees:	031002060302
# **** Recommend changing COMID to 16918094
# 
# 
# ===***===***===***===***===***===***===***===***===***===***===***===***===***===***
#   coordinates  site_no    comid        huc12 decade       lon      lat    n
# (558732.2, 1350551) 07032200 14199483 080102110102   2000 -89.81898 35.04981 3653
# 07032200 NONCONNAH CREEK NEAR GERMANTOWN, TN
# The assigned COMID 14199483 does not exist.  Looking through the COMID list,
# RRK thinks this was a typo, albeit incorrect for the site location.
# COMID 14199473 is on the same river, but also not correct.
# Gage is located on COMID: 14199477
# Most DS COMID in HUC12: 14199375
# HUC12 agrees: 080102110102
# **** Recommend changing COMID to 14199477
# 
# 
# ===***===***===***===***===***===***===***===***===***===***===***===***===***===***
#   coordinates  site_no    comid        huc12 decade       lon      lat    n
# (535831.8, 1482204) 07045000   766856 080202040615   1960 -89.97949 36.23652 3653
# 07045000 LITTLE RIVER DITCH 66 NEAR KENNETT MO
# This is a ditch, and the gage plots between four parallel lines (ditches). There is
# another Little River Kennett, MO (see section of '1/3' log10 cycle at end of this file.)
# Nearest COMID: 766886 (Little River)
# COMID at DS HUC12: 	766886 (same)
# HUC12 agrees: 080202040615 (long skinny HUC12 between other ditches)
# Original COMID is 2 ditches to the east
# **** Recommend changing COMID to 766886
# 
# 
# ===***===***===***===***===***===***===***===***===***===***===***===***===***===***
#   coordinates  site_no    comid        huc12 decade       lon      lat    n
# (505060.2, 1398264) 07047600  3672364 080202031210   1960 -90.38010 35.50508 3653
# 07047600 TYRONZA RIVER NEAR TYRONZA, ARK
# Gage is located on COMID: 3672364
# Most DS COMID for HUC12: 3674276
# This gage is correct --- Tentative conclusion is that NWIS then is erroneous.
# HUC12 agrees: 080202031210
# **** Recommend no change to COMID.
# 
# 
# ===***===***===***===***===***===***===***===***===***===***===***===***===***===***
#   coordinates  site_no    comid        huc12 decade       lon      lat    n
# (525437.4, 1215141) 07281000 18046596 080302020600   1970 -90.27648 33.85984 3652
# 07281000 TALLAHATCHIE RIVER AT SWAN LAKE, MS
# wonky hydrology
# Different HUC12: HUC12	080302020602
# nearest COMID of same name: ComID	18050376
# most DS COMID for HUC12: 18047158
# **** Recommend changing COMID to 18050376
# 
# 
# ===***===***===***===***===***===***===***===***===***===***===***===***===***===***
#   coordinates  site_no    comid        huc12 decade      lon      lat    n
# (392675.5, 1054107) 07368000 19350501 080500010801   1950 -91.7979 32.48125 3652
# 07368000 Boeuf River near Girard, LA
# Gage is located on COMID 19350501
# Most DS COMID forthe HUC12: 19350501
# This gage is correct --- Tentative conclusion is that NWIS then is erroneous. This
# site does not capture the entire basin because of channels intertwining as a result
# of low-gradient hydrology.
# HUC12 agrees: 080500010801
# **** Recommend no change to COMID
# 
# 
# ===***===***===***===***===***===***===***===***===***===***===***===***===***===***
#   coordinates  site_no    comid        huc12 decade       lon  lat    n
# (381364.9, 1055706) 07369000 19350633 080500011304   1950 -91.91806 32.5 3652
# 07369000 Bayou Lafourche near Crew Lake, LA
# site located on COMID: 19350633
# DS HUC12 COMID:  19350645
# This gage is correct --- Tentative conclusion is that NWIS then is erroneous. This
# site does not capture the entire basin because of channels intertwining as a result
# of low-gradient hydrology.
# HUC12	agrees: 080500011304
# **** Recommend no change to COMID
# 
# 
# ===***===***===***===***===***===***===***===***===***===***===***===***===***===***
#   coordinates  site_no   comid        huc12 decade       lon     lat    n
# (41828.68, 1198606) 07335400 8349005 111401010605   1970 -95.54468 33.8526 3652
# Sanders Creek near Chicota, Texas
# Gage appears installed at about the same time of reservoir construction. Suspect
# # differences in treatment of segments and inserted "artificial paths" somehow is
# # to blame. This site could be a real test for NHDPlus (or our understanding) whether
# # they in fact have a COMID 8348979. By "real test," WHA means that we are addressing
# # things the same.
# HUC12	agrees: 111401010605
# Nearest COMID: COMID	8348979 ("artificial path," first segment above but touching the dam)
# Most DS COMID for HUC12: 8348979
# **** Recommend changing COMID to 8348979
# 
# 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   +   NEW SECTION. These are outliers of |>1/3| log10 and |<1/2| log10.
# +   We had five others but they were all small watershed and are not further considered
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   ===***===***===***===***===***===***===***===***===***===***===***===***===***===***
#   coordinates  site_no    comid        huc12 decade       lon      lat    n
# (1135598, 870408) 02326900 10323850 031200011002   1960 -84.14825 30.27143 3653
# 02326900 ST. MARKS RIVER NEAR NEWPORT, FLA. (535 sq mi)
# nearest COMID: 10323850
# DS COMID for HUC12: 10320308
# The NHD line work crosses back and forth over the same HUC boundary several times.
# This is the COMID for the first crossing)
# HUC12	031200011003
# **** Recommend changing COMID to 10323850
# 
# 
# ===***===***===***===***===***===***===***===***===***===***===***===***===***===***
#   coordinates  site_no    comid        huc12 decade       lon      lat    n
# (1388758, 691638.1) 02310800 16954130 031002080102   1960 -81.81869 28.36084 3653
# 02310800 WITHLACOOCHEE RIVER NR EVA, FLA. (130 sq. mi)
# nearest COMID: 16954130
# DS COMID for HUC12:   16954130
# HUC12	031002080102  (site sits right on HUC12 boundary - reporting the upstream HUC12)
# **** Recommend changing COMID to 16954130
# 
# 
# ===***===***===***===***===***===***===***===***===***===***===***===***===***===***
#   coordinates  site_no    comid        huc12 decade       lon      lat    n
# (535931.7, 1482190) 07044000 766888 080202040615   1950 -89.97838 36.23633 3652
# 07044000 Little River Ditch 251 near Kennett, MO (883 sq mi)
# 1 of 4 parallel ditches (similar to site at top of file), long skinny HUC12
# nearest COMID: 766886
# DS COMID for HUC12:   767886
# HUC12	080202040615
# **** Recommend changing COMID to 766886