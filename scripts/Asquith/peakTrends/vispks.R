
# Must go get this non-CRAN package.
# https://github.com/wasquith-usgs/MGBT
# To install and run in R...
# install.packages("devtools")
# devtools::install_github("wasquith-usgs/MGBT")
# install.packages(c("dataRetrieval", "lmomco")) # These are on the CRAN.

library(lmomco) # though MGBT will load it
library(MGBT)
library(dataRetrieval)
library(feather)

pdffile <- "vispk.pdf"
rdfile  <- "vispk.RData"
sdecade <- 1890
edecade <- 2010

# If these are defined once, then one can comment out the next two lines
# so that development or testing on (say) plotPeaks() could be made w/o
# repulling the data from the Internet because we have a caching mechanism.
PK  <- new.env() # environment to house ALL peaks
LOT <- new.env() # environment to house the low-outlier thresholds
# One accessed the information in the environment with syntax like:
# get.that.data.frame <- PK$"08167000"
# The LOT is cached too as the MGBT() process can be time consuming.
sites <- NULL

if(! is.na(rdfile) & file.exists(rdfile)) {
  load(rdfile) # We can have a premade RData file laying around instead
  sites <- ls(PK)
} else if(! is.na(rdfile)) { # decade1950plus_site_list.csv
  D <- read_feather("../../../data/gage/all_gage_data.feather")
  sites <- unique(D$site_no)
} else {
  stop("No data available for processing")
}


SF <- dataRetrieval::readNWISsite(sites) # grab the site file
SF <- SF[SF$agency_cd == "USGS",] # only retain USGS, Texas has a couple of
# gages cross listed as USCOE, so we need to drop potential duplicates.
if(length(sites) != length(SF$station_nm)) {
  stop("fatal mismatch in site number and station name lengths")
}
idname <- paste(sites, SF$station_nm) # for a titling of the plots


n <- length(sites); i <- 0 # a counter
allpv <- alltau <- sompv <- somtau <- n.tau <- n.nonlot <- missing <- rep(NA, n)
if(! is.na(pdffile)) pdf(pdffile, useDingbats=FALSE)
for(site in sites) {
  i <- i + 1
  message("  site=",site," as ",i," of ",n)
  pk <- lot <- NULL # insurance policy to NULLify if code become "complicated"
  if(exists(site, PK)) { # trigger the cache if present
    pk <- get(site, envir=PK) # by this test, the user could avoid repulling
    if(is.null(pk)) {
      message(" missing ", site)
      missing[i] <- site
    }
  } else {
    pk <- dataRetrieval::readNWISpeak(site,convert=FALSE) # CONVERT==FALSE!
    if(is.null(pk) | length(pk) == 1) {
      assign(site, NULL, envir=PK) # store the peak tables separately
      message(" missing ", site)
    } else {
      pk$peak_va <- as.numeric(pk$peak_va) # because convert=FALSE
      pk <- MGBT::splitPeakCodes(MGBT::makeWaterYear(pk)) # see documentation
      pk$idname <- idname[sites == site]
      # these operations are added convenience columns to the table retrieved
      assign(site, pk, envir=PK) # store the peak tables separately
    }
  }
  if(exists(site, LOT)) { # trigger the cache if present
    lot <- get(site, envir=LOT) # again a caching mechanism
  } else {
    if(is.null(pk) | length(pk) == 1) {
      assign(site, 0, envir=LOT)
    } else {
      lot <- MGBT::MGBT(pk$peak_va[! is.na(pk$peak_va)])$LOThresh
      assign(site, lot, envir=LOT)
    }
  }
  if(is.null(pk) | length(pk) == 1) {

  } else {
    plotPeaks(pk, lot=lot, codes=TRUE,
                  ylab="Peak discharge, in cubic feet per second",
                  xlim=c(sdecade,edecade), site=pk$idname[1])
    pdf(paste0(site,".pdf"), useDingbats=FALSE)
      plotPeaks(pk, lot=lot, codes=TRUE,
                    ylab="Peak discharge, in cubic feet per second",
                    xlim=c(sdecade,edecade), site=pk$idname[1])
    dev.off()
    tmpa <- pk[pk$appearsSystematic,]
    opts <- options(warn=-1)
    all  <- cor.test(tmpa$water_yr, tmpa$peak_va, method="kendall")
    alltau[i] <- all$estimate; allpv[i] <- all$p.value; n.tau[i] <- length(tmpa$water_yr)
    som <- NULL
    if(lot != 0) {
      tmpb <- tmpa[tmpa$peak_va > lot, ]
      som <- cor.test(tmpb$water_yr, tmpb$peak_va, method="kendall")
      somtau[i] <- som$estimate; sompv[i] <- som$p.value; n.nonlot[i] <- length(tmpb$water_yr)
    }
    options(opts)
    # If you desire to test or change plotPeaks(), open its sources, changes the
    # function name to plotPeaks2 and change it here in this loop, then modify
    # and source it, and rerun this script. Because of this potential desire,
    # the pkgcoloncolon is not prepended. These are prepended to immediately show which
    # package is provide what function even though this strictly is not needed.
  }
}
if(! is.na(pdffile)) dev.off()


TAU <- data.frame(site_no=sites, n.tau= n.tau, tau=alltau, tau_p.value=allpv,
                  n.nonlot=n.nonlot, nonlot_tau=somtau, nonlot_p.value=sompv,
                  sign_reinforcement=as.factor(sign(alltau) + sign(somtau)))

missing <- missing[! is.na(missing)]

if(! is.na(rdfile)) save(PK, LOT, SF, TAU, missing, file=rdfile)

library(akqdecay)
DV <- new.env()
fill_dvenv(siteNumbers=missing, envir=DV)
DVMX <- NULL
for(site in sort(ls(DV))) {
  tmp <- get(site, envir=DV)
  h <- aggregate(tmp, by=list(tmp$wyear), max)
  n <- aggregate(tmp, by=list(tmp$wyear), function(t) length(t))
  d <- data.frame(agency_cd=h$agency_cd, site_no=h$site_no, peak_dt=rep("DATE",length(h$site_no)),
                  peak_tm=NA, peak_va=h$Flow, peak_cd="1", gage_ht=NA, gage_ht_cd=NA, year_last_pk=NA,
                  ag_dt=NA, ag_tm=NA, ag_gage_ht=NA, ag_gage_ht_cd=NA,
                  water_yr=h$wyear, num_dvs=n$wyear, stringsAsFactors=FALSE)
  d <- d[d$num_dvs >= 357, ] # reject year if about 7 days missing
  if(is.null(DVMX)) {
    DVMX <- d
  } else {
    DVMX <- rbind(DVMX,d)
  }
}
DVMX$peak_cd[DVMX$site_no == "02301802"] <- "1,C" # 02301802 TAMPA BYPASS CANAL AT S−160,AT TAMPA FL
DVMX$peak_cd[DVMX$site_no == "02304500"] <- "1,6,C" # 02304500 HILLSBOROUGH RIVER NEAR TAMPA FL
DVMX$peak_cd[DVMX$site_no == "02343500"] <- "1,C" # 02343500 CHATTAHOOCHEE R AT COLUMBIA, AL
DVMX$peak_cd[DVMX$site_no == "07024900"] <- "1" # 07024900 RUTHERFORD FORK OBION RIVER NEAR MILAN, TENN
DVMX$peak_cd[DVMX$site_no == "07359001"] <- "1,6" # 07359001 Ouachita River below Remmel Dam at Jones Mill, AR

for(i in 1:length(DVMX$peak_va)) {
  site <- DVMX$site_no[i]; wyear <- DVMX$water_yr[i]; peak  <- DVMX$peak_va[i]
  tmp <- get(site, envir=DV); tmp <- tmp[tmp$wyear == wyear,]
  dt <- as.character(tmp$Date[tmp$Flow == peak])
  if(length(dt) == 0) {
    message(site, " ", peak, " ", wyear)
  }
  DVMX$peak_dt[i] <- dt[length(dt)]
}
write.table(DVMX, file="noNWISpeaks-sitesIn-all_gage_data_feather.txt", quote=FALSE, sep="\t", row.names=FALSE)
for(site in sort(unique(DVMX$site_no))) {
  tmp <- DVMX[DVMX$site_no == site,]
  tmp$peak_va <- as.numeric(tmp$peak_va) # because convert=FALSE
  tmp <- MGBT::splitPeakCodes(MGBT::makeWaterYear(tmp)) # see documentation
  tmp$idname <- idname[grep(site,idname)]
  assign(site, tmp, envir=PK)
}

#plot(TAU$tau, TAU$nonlot_tau, col=abs(TAU$nonlot_p.value < 0.05)+abs(TAU$tau_p.value < 0.05)+2)
