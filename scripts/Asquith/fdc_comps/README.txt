
This README.txt file accompanying this data release demonstrates in detail
the core algorithm to compute the statistics reported. The audience intended
are those developers with considerable understanding of the idioms of the
R statistical programming language as well as the installation of add-on
packages from the Comprehensive R Archive Network (CRAN). The demonstration
is expected to run "out of the box," but it does not provide the facilities
for large-scale looping through a list of sites and attendant error trapping.
Reading this README.txt with a text editor and a monospaced font is suggested
to enhance the organization of the code, comments, and results.


The demonstrated components of the algorithm compute the decadal flow-
duration curve (FDC) quantiles, the minimum and maximum of the flows (zeros
included), the L-moments of the nonzero portion of the FDC, the median of the
nonzero portion, and the fraction of decadal zero flow. The demonstration also
shows the treatment for missing record. The fdclmr_decade() defined below also
provides for a handling of logarithmic transformation, which was not used for
this data release.


The demonstration uses a selected U.S. Geological streamflow-gaging station
(streamgage) and subsequent trimming for the decadal range 1950--2000, and is
at the end of this README.txt file. The streamgage is included in this data
release. The streamgage was chosen because it has some decades with zeroflow
values and some without so that the breadth of the algorithm is demonstrated.
An Internet connection is required because the demonstration dynamically will
pull the period of record of daily-mean streamflow values for the streamgage.


# -----------------------------------------------------------------------------
#  Basic Algorithm for U.S. Geological Streamflow-Gaging Station 08167000.
#  The text that follows can be cut and pasted into an R language console and
#  is expected to function if two packages are already installed. Selected
#  (for brevity) expected results are shown in the last lines of this file. 
# -----------------------------------------------------------------------------
"fdclmr_decade" <-
function(dvtable, missing.days=7, site="", log=FALSE, subzero=NULL, plusit=1) {
  if(length(unique(dvtable$site_no)) > 1) {
    warning("can not have more than one streamgage in the daily value table, ",
            "returning NA immediately"); return(NA)
  }
  site <- as.character(site[1]) # site is a special override on the site id

  if(length(dvtable$year) == 0) return(NA)

  # selected flow-duration curve probabilities
  probs <- c(0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 25, 30, 40, 50, 60,
             70, 75, 80, 90, 95, 98, 99, 99.5, 99.8, 99.9, 99.95, 99.98) / 100
  dds <- unique(dvtable$decade)
  avail <- 365*10 - missing.days*10 # missing record allowance
  nlmoms <- 8 # number of L-moments
  empty <- list(fdc_quantiles=rep(NA, length(probs)),  min=NA, max=NA, pplo=NA,
                lambdas=rep(NA,nlmoms),ratios=rep(NA,nlmoms), median_nonzero=NA,
                source="not valid")
  zz <- data.frame(site_no=site, decade=NA, n=NA, nzero=NA, pplo=NA,
          min=NA,  f0.02=NA, f0.05=NA, f0.1=NA, f0.2=NA, f0.5=NA, f01=NA,
          f02=NA, f05=NA, f10=NA, f20=NA, f25=NA, f30=NA, f40=NA, f50=NA,
          f60=NA, f70=NA, f75=NA, f80=NA, f90=NA, f95=NA, f98=NA, f99=NA,
          f99.5=NA, f99.8=NA, f99.9=NA, f99.95=NA, f99.98=NA, max=NA,
          L1=NA, L2=NA, T3=NA, T4=NA, T5=NA, T6=NA, T7=NA, T8=NA,
          median_nonzero=NA); zz <- zz[0,]

  for(d in dds) {
     fdc <- dvtable$Flow[dvtable$decade == d] # grab flow for a decade
     nzero <- length(fdc[fdc == 0])
     n <- length(fdc[! is.na(fdc)]); fdc <- fdc[! is.na(fdc)]
     if(log) { # if a user desires logarithms (experimental)
        nzero <-        length(fdc[fdc == 0])
        if(! is.null(subzero)) fdc[fdc == 0] <- subzero
        if(! is.null(plusit )) fdc <- fdc + plusit
        opts <- options(warn=-1)
          fdc <-            log10(fdc)
          if(length(fdc[is.nan(   fdc)]) > 0) message("NaN ",site," for ", d)
          fdc <-    fdc[is.finite(fdc)]
        options(opts)
     }
     if(any(is.na(fdc))) { # next few lines deal with some error trapping
        message("a least one missing value for year or decade ",d)
        lmr <- empty
     } else if (length(unique(fdc)) == 1) {
        lmr <- empty
     } else {
        if(n <= avail) {
           lmr <- empty
        } else {
           # See the documentation of the lmomco package for more information.
           fdclo <- lmomco::x2xlo(fdc) # L-moments of zero left-truncation
           #print(fdclo$xin) # The xin are the values "left in" or nonzero.
           lmr <-         lmomco::lmoms(fdclo$xin, nmom=nlmoms, no.stop=TRUE)
           lmr$median_nonzero <- median(fdclo$xin) # median of the nonzero part
           if(! lmomco::are.lmom.valid(lmr)) { # bailout to probability
              # weighted moment by plottingposition L-moments.
              # This ensures maximal availability of the L-moments.
              lmr <- lmomco::pwm2lmom(lmomco::pwm.pp(fdclo$xin, nmom=8))
           }
           lmr$pplo <- fdclo$pp # percentage of zero flow
           lmr$min <- min(fdc); lmr$max <- max(fdc) # minimum and maximums
           lmr$fdc_quantiles <- quantile(fdc, probs=probs, type=6)
        }
     }
     q <- lmr$fdc_quantiles
     tmp <- data.frame(site_no=site, decade=d, n=n, nzero=nzero, pplo=lmr$pplo,
              min=lmr$min, f0.02=q[1], f0.05=q[2], f0.1=q[3], f0.2=q[4],
              f0.5=q[5], f01=q[6],  f02=q[7],  f05=q[8],  f10=q[9],
              f20=q[10], f25=q[11], f30=q[12], f40=q[13], f50=q[14],
              f60=q[15], f70=q[16], f75=q[17], f80=q[18], f90=q[19],
              f95=q[20], f98=q[21], f99=q[22], f99.5=q[23], f99.8=q[24],
              f99.9=q[25], f99.95=q[26], f99.98=q[27], max=lmr$max,
              L1=lmr$lambdas[1], L2=lmr$lambdas[2], T3=lmr$ratios[3],
              T4=lmr$ratios[4],  T5=lmr$ratios[5],  T6=lmr$ratios[6],
              T7=lmr$ratios[7],  T8=lmr$ratios[8],
              median_nonzero=lmr$median_nonzero, stringsAsFactors=FALSE)
     zz <- rbind(zz, tmp); row.names(zz) <- NULL
  }
  return(zz) # one to many rowed data.frame depending decades available
}

library(dataRetrieval); library(lmomco)

siteNumber <- "08167000" # selected streamgage for demonstration

# See the documentation of the dataRetrieval package for more information.
dvs <- dataRetrieval::readNWISdv(siteNumber, startDate="", endDate="",
                                 parameterCd="00060", statCd="00003")
dvs <- dataRetrieval::renameNWISColumns(dvs)

dvs$year <- sapply(strsplit(as.character(dvs$Date), split="-"),
                   function(i) return(i[1])) # extract the year
dvs$decade <- as.integer(as.numeric(sub("\\d$", "", dvs$year)) * 10)

fdc <- fdclmr_decade(dvs, site=dvs$site_no[1]) # compute statistics

# The subselection for >= 1950 was used for this data release. The secondary
# subselection for <= 2000 was not needed because at the time of release, the
# 2010 decade was incomplete.
#                    (>=, greater than or equal; <=, less than or equal)
fdc <- fdc[fdc$decade >= 1950    &    fdc$decade <= 2000,]

head(fdc[,1:10]) # inspect the results as shown below for the first ten columns.
# 
#   site_no decade    n nzero       pplo  min f0.02    f0.05    f0.1   f0.2
#  08167000   1950 3652   278 0.07610183  0.0   0.0  0.00000  0.0000  0.000
#  08167000   1960 3653    48 0.01313629  0.0   0.0  0.00000  0.0000  0.000
#  08167000   1970 3652     0 0.00000000  5.8   5.8  5.80000  6.0000  8.100
#  08167000   1980 3653     0 0.00000000  1.8   1.8  1.80000  2.3270  4.554
#  08167000   1990 3652     0 0.00000000  9.7   9.7 10.77445 11.0000 13.306
#  08167000   2000 3653     0 0.00000000 12.0  12.0 12.16540 12.6616 13.100
#
# The results show the decadal fraction of zeroflow as 0.013 for the 1960
# decade and the minimum daily-mean streamflow as 12 cubic feet per second
# for 2000 decade.
