A brief description of the master lists of site identification numbers for U.S. Geological Survey streamflow-gaging stations. The lists referred to here are exclusive to only the site number. Tables with data, watershed properties, landuse, and other metadata for purposes of analysis are listed elsewhere in this archive.

file =  full_site_list.csv
These are the 1,379 streamgages used to make the query of all daily values of streamflow for the period of record through 2016-09-30.
__END__


file = decade1950plus_site_list.csv
These are the 1,033 streamgages having "complete" decades (10 * one week = 70 days permitted to be missing) for decades after and including 1950. Exploratory analyses for this project indicates that core statistical regionalization needs to be restricted to at least the "modern" period after World War II. This criteria of decade and post 1950 decided on 12/8/2017 by SCW, WHA, and RRK. The FDC data for these can be found in the file fdc_lmr_pplo.feather.
__END__

file = decade1950plus_plus2010-15_site_list.csv
These are the 1,166 streamgages within the union of sites between files fdc_lmr_pplo.feather and fdc_lmr_pplo2010-15.feather. Special accomodation for a 2010 partial decade was made for the 6 years between 2010-15. The total 1,166 is 1,033 + 133. There are a substantial number of very young gages in the study area that uniquely showup only in the 2010 partial decade (about 80) so the 133 is a substantial fraction of the 1,033 sites for the decades 1950--2000.

> library(feather)
> setwd("~/Projects/RESTORE/restore-2018/data/gage")
> fdc10 <- read_feather("fdc_lmr_pplo2010-15.feather")
> fdc00 <- read_feather("fdc_lmr_pplo.feather")
> sites <- sort(unique(c(fdc00$site, fdc10$site)))
> length(sites)
[1] 1166
> site_list <- sort(unique(c(fdc00$site, fdc10$site)))
> site_list <- data.frame(site_no=site_list, stringsAsFactors =FALSE)
write.table(site_list, file="decade1950plus_plus2010-15_site_list.csv", quote=TRUE, row.names=FALSE)

__END__


