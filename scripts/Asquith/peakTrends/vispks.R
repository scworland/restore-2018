
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
allpv <- alltau <- sompv <- somtau <- n.tau <- n.nonlot <- rep(NA, n)
if(! is.na(pdffile)) pdf(pdffile, useDingbats=FALSE)
for(site in sites) {
  i <- i + 1
  message("  site=",site," as ",i," of ",n)
  pk <- lot <- NULL # insurance policy to NULLify if code become "complicated"
  if(exists(site, PK)) { # trigger the cache if present
    pk <- get(site, envir=PK) # by this test, the user could avoid repulling
    if(is.null(pk)) message(" missing ", site)
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
                  n.nonlot=n.nonlot, nonlot_tau=somtau, nonlot_p.value=sompv)

if(! is.na(rdfile)) save(PK, LOT, SF, TAU, file=rdfile)


#plot(TAU$tau, TAU$nonlot_tau, col=abs(TAU$nonlot_p.value < 0.05)+abs(TAU$tau_p.value < 0.05)+2)
