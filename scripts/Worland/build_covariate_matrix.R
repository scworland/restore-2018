
library(tidyverse)
library(lubridate)
library(feather)
library(geoknife)

summarize = dplyr::summarise
query = geoknife::query

# query available variables in prism dataset
daymet = webdata(url='https://thredds.daac.ornl.gov/thredds/dodsC/daymet-v3-agg/na.ncml')
query(daymet,'variables')

# load data
precip <- webdata(list(
  times = as.POSIXct(c('1970-01-01','2010-12-01')),
  url = 'https://thredds.daac.ornl.gov/thredds/dodsC/daymet-v3-agg/na.ncml',
  variables = 'prcp'))

tmax <- webdata(list(
  times = as.POSIXct(c('1970-01-01','2010-12-01')),
  url = 'https://thredds.daac.ornl.gov/thredds/dodsC/daymet-v3-agg/na.ncml',
  variables = 'tmax'))

tmin <- webdata(list(
  times = as.POSIXct(c('1970-01-01','2010-12-01')),
  url = 'https://thredds.daac.ornl.gov/thredds/dodsC/daymet-v3-agg/na.ncml',
  variables = 'tmin'))

vp <- webdata(list(
  times = as.POSIXct(c('1970-01-01','2010-12-01')),
  url = 'https://thredds.daac.ornl.gov/thredds/dodsC/daymet-v3-agg/na.ncml',
  variables = 'vp'))

swr <- webdata(list(
  times = as.POSIXct(c('1970-01-01','2010-12-01')),
  url = 'https://thredds.daac.ornl.gov/thredds/dodsC/daymet-v3-agg/na.ncml',
  variables = 'srad'))

# create huc12 stencil


stencil <- webgeom(geom='upload:huc12alb', attribute="HUC_12")
HUCs <- query(stencil, 'values')

query(webgeom(geom='upload:huc12alb'), 'geoms')

# download data
precip.job <- geoknife(stencil=all.cnty, fabric=precip, wait = F, email = 'scworland@usgs.gov')
tmax.job <- geoknife(stencil=all.cnty, fabric=tmax, wait = F, email = 'scworland@usgs.gov')
tmin.job <- geoknife(stencil=all.cnty, fabric=tmin, wait = F, email = 'scworland@usgs.gov')




