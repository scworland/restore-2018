library(feather)
library(proxy)
library(akqdecay)
library(mgcv)
library(sp)
library(rgeos)
library(lmomco)
library(Lmoments)
library(survival)
library(cenGAM)
library(hydroGOF)

source("gamIntervals.R")


load("../../../../GIS/RESTORE_MGCV_BND.RData") # "RESTORE_MGCV_BND.RData"

FDC <- read_feather("../../../data/gage/all_gage_flow_stats.feather")
#nm <- names(FDC); nm[5] <- "dec_long_va"; nm[6] <- "dec_lat_va"; names(FDC) <- nm
#write_feather(FDC, "../../../data/gage/all_gage_flow_stats.feather")

COV <- read_feather("../../../data/gage/all_gage_covariates.feather")
#nm <- names(COV); nm[4] <- "dec_long_va"; nm[5] <- "dec_lat_va"; names(COV) <- nm
#write_feather(COV, "../../../data/gage/all_gage_covariates.feather")

SO <- read_feather("../../../data/gage/all_gage_solar.feather")
COV <- cbind(COV, SO)

FDC$key <- paste(FDC$site_no,":",FDC$decade, sep="")
COV$key <- paste(COV$site_no,":",COV$decade, sep="")
suppressWarnings(DD <- merge(FDC, COV, by="key"))
DD$key <- NULL
#Warning message:
#In merge.data.frame(FDC, COV, by = "key") :
#  column names ‘comid.y’, ‘site_no.y’, ‘huc12.y’, ‘dec_long_va.y’, ‘dec_lat_va.y’, ‘decade.y’ are duplicated in the result
sum(DD$comid.x       != DD$comid.y      )
sum(DD$site_no.x     != DD$site_no.y    )
sum(DD$huc12.x       != DD$huc12.y      )
sum(DD$dec_long_va.x != DD$dec_long_va.y)
sum(DD$dec_long_va.x != DD$dec_long_va.y)
sum(DD$decade.x      != DD$decade.y     )
nm <- names(DD); nm[1] <- "comid"; nm[2] <- "site_no"; nm[3] <- "huc12"; nm[4] <- "decade"
nm[5] <- "dec_long_va"; nm[6] <- "dec_lat_va"
names(DD) <- nm
DD$comid.y <- NULL; DD$site_no.y <- NULL; DD$huc12.y <- NULL
DD$dec_long_va.y <- NULL; DD$dec_lat_va.y <- NULL; DD$decade.y <- NULL
length(DD$site_no); length(unique(DD$site_no))

DD[! complete.cases(DD),] # Better be zero rows!


sites <- unique(DD$site_no)
sitefile <- dataRetrieval::readNWISsite(sites)
sitefile <- sitefile[sitefile$agency_cd != "USCE",]

CDA <- sitefile$contrib_drain_area_va
CDA[is.na(CDA)] <- sitefile$drain_area_va[is.na(CDA)]
CDA <- pmin(sitefile$drain_area_va, sitefile$contrib_drain_area_va, na.rm=TRUE)
sitefile$contrib_drain_area_va <- CDA <- CDA*2.589988 # move to km2
#write.table(sitefile, file="sitefile.txt", quote=FALSE, row.names=FALSE, sep="\t")

LATLONG <- paste0("+proj=longlat +ellps=GRS80 ",
                  "+datum=NAD83 +no_defs +towgs84=0,0,0")
LATLONG <- sp::CRS(LATLONG)
ALBEA <- paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 ",
                "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
ALBEA <- sp::CRS(ALBEA)

SF <- data.frame(site_no=sitefile$site_no,
                 CDA=log10(CDA), stringsAsFactors=FALSE)
DD <- merge(DD, SF, all=TRUE)


DD <- SpatialPointsDataFrame(cbind(DD$dec_long_va, DD$dec_lat_va), DD,
                             proj4string=LATLONG)
DD <- spTransform(DD, ALBEA)
XY <- coordinates(DD); DD$east <- XY[,1]/1000; DD$north <- XY[,2]/1000; rm(XY)
DD$x <- DD$east; DD$y <- DD$north

ix <- length(2:(length(bnd_poly_aea[,1])-1))
ix <- c(1,sort(sample(ix, size=20000, replace=FALSE)),length(bnd_poly_aea[,1]))
bnd <- list(x=bnd_poly_aea[,1]/1000, y=bnd_poly_aea[,2]/1000)
bnd <- list(bnd)
bndthin <- list(x=bnd_poly_aea[ix,1]/1000, y=bnd_poly_aea[ix,2]/1000)
bndthin <- list(bndthin)
bnd <- bndthin
plot(bnd[[1]]$x,bnd[[1]]$y,type="l", col=8, lwd=.6)
points(DD$east, DD$north, pch=4, lwd=.5, cex=0.9, col=rgb(0,.5,0.5,.2))

# This could be one way to build a grid of knots. However, a
# large number do not seem to be necessary. Further, knots in
# areas with limited data cause poor results.
w <- seq(-1000, 1500, by=200)
v <- seq(  300, 1700, by=200)
m <- length(w); n <- length(v)
w <- rep(w, n)
h <- NULL
for(i in 1:n) { h <- rbind(h,rep(v[i], m)) }
x <- w; y <- as.vector(h)
ind <- mgcv::inSide(bnd,x,y)
knots <- data.frame(x=x[ind], y=y[ind])
rm(x,y,v,w,h)
points(knots$x, knots$y, pch=16, cex=0.3, col=2)
text(knots$x, knots$y, row.names(knots))
# Another problem with this approach is that inevitiably we
# get knots inside the boundary according to inSide() but when
# it comes time to actually build the soapfilm smoother, we get
# errors of knots outside the boundary when they seem inside.
# Internally, the mgcv logic must have some type of fuzzy test but
# I have not found any documentation.
knots <- knots[c(20, 24),]
#knots <- knots[-c(1, 2, 4, 8, 9, 12, 13, 21),]
x <- knots$x; y <- knots$y
#knots <- data.frame(x=c(x, -200), y=c(y, 800)); rm(x,y)
points(knots$x, knots$y, pch=16, cex=1.1, col=4)

length(DD$site_no)

DD$decade       <- as.factor(DD$decade);       levels(DD$decade)
DD$cat_soller   <- as.factor(DD$cat_soller);   levels(DD$cat_soller)
DD$soller       <- as.factor(DD$soller);       levels(DD$soller)
DD$cat_aquifers <- as.factor(DD$cat_aquifers); levels(DD$cat_aquifers)
DD$aquifers     <- as.factor(DD$aquifers);     levels(DD$aquifers)
DD$bedperm      <- as.factor(DD$bedperm);      levels(DD$bedperm)
DD$cat_physio   <- as.factor(DD$cat_physio);   levels(DD$cat_physio)
DD$physio       <- as.factor(DD$physio);       levels(DD$physio)
DD$cat_ecol3    <- as.factor(DD$cat_ecol3);    levels(DD$cat_ecol3)
DD$ecol3        <- as.factor(DD$ecol3);        levels(DD$ecol3)
DD$hlr          <- as.factor(DD$hlr);          levels(DD$hlr)
DD$statsgo      <- as.factor(DD$statsgo);      levels(DD$statsgo)
DD$ed_rch_zone  <- as.factor(DD$ed_rch_zone);  levels(DD$ed_rch_zone)
unique(DD$site_no[DD$ed_rch_zone == "1"])

# only nodata for cat_soller, cat_aquifers, and cat_physio have 1 nodata : site_no 02359315
summary(DD$cat_soller)
summary(DD$soller)
summary(DD$cat_aquifers)
summary(DD$aquifers)
summary(DD$bedperm)
summary(DD$cat_physio)
summary(DD$physio)
summary(DD$cat_ecol3)
summary(DD$ecol3)
summary(DD$hlr)
summary(DD$statsgo)
summary(DD$ed_rch_zone)



DD$ppt_mean      <- log10(DD$ppt_mean)
DD$temp_mean     <- log10(DD$temp_mean)
DD$basin_area    <- log10(DD$basin_area)
DD$basin_slope   <- log10(DD$basin_slope/100)
flood_storage_offset <- 1E-6; # first even log10 cycle below min(DD$flood_storage[DD$flood_storage > 0])
DD$flood_storage <- log10(DD$flood_storage + flood_storage_offset)

# Transformation and Retransformation Functions for the Sin Transformation of
# Percentile data
dotransin <- function(p) 2*asin(sqrt(p/100))
retransin <- function(p)    sin(p/2)^2*100

DD$barren              <-  dotransin(DD$barren)
DD$cultivated_cropland <-  dotransin(DD$cultivated_cropland)
DD$deciduous_forest    <-  dotransin(DD$deciduous_forest)
DD$developed           <-  dotransin(DD$developed)
DD$evergreen_forest    <-  dotransin(DD$evergreen_forest)
DD$grassland           <-  dotransin(DD$grassland)
DD$hay_pasture         <-  dotransin(DD$hay_pasture)
DD$herbaceous_wetland  <-  dotransin(DD$herbaceous_wetland)
DD$mixed_forest        <-  dotransin(DD$mixed_forest)
DD$shrubland           <-  dotransin(DD$shrubland)
DD$water               <-  dotransin(DD$water)
DD$woody_wetland       <-  dotransin(DD$woody_wetland)
DD$bfi                 <-  dotransin(DD$bfi)

#DD$edwards_rechzone <- 0
#DD$edwards_rechzone[DD$site_no == "08155300"] <- 1
#DD$edwards_rechzone[DD$site_no == "08155400"] <- 1
#DD$edwards_rechzone[DD$site_no == "08156800"] <- 1
#DD$edwards_rechzone[DD$site_no == "08181400"] <- 1
#DD$edwards_rechzone[DD$site_no == "08184000"] <- 1
#DD$edwards_rechzone[DD$site_no == "08185000"] <- 1
#DD$edwards_rechzone[DD$site_no == "08190500"] <- 1
#DD$edwards_rechzone[DD$site_no == "08197500"] <- 1
#DD$edwards_rechzone[DD$site_no == "08198500"] <- 1
#DD$edwards_rechzone[DD$site_no == "08200700"] <- 1
#DD$edwards_rechzone[DD$site_no == "08202700"] <- 1
#DD$edwards_rechzone <- as.factor(DD$edwards_rechzone)

DDo <- DD

pdf("first_diagnostic_on_area.pdf", useDingbats=FALSE, width=7.5, height=6.5)
par(mgp=c(3,0.5,0)) # going to tick to the inside, change some parameters
plot(10^DDo$CDA, 10^DDo$basin_area, lwd=0.5, xaxt="n", yaxt="n", log="xy",
     xlab="", ylab="", xlim=c(1,300000), ylim=c(1,300000))
add.log.axis(side=2,      tcl=0.8*abs(par()$tcl), two.sided=TRUE)
add.log.axis(logs=c(1),   tcl=-0.5*abs(par()$tcl), side=2, two.sided=TRUE)
add.log.axis(logs=c(1),   tcl=+1.3*abs(par()$tcl), side=2, two.sided=TRUE)
add.log.axis(side=1,      tcl=0.8*abs(par()$tcl), two.sided=TRUE)
add.log.axis(logs=c(1),   tcl=-0.5*abs(par()$tcl), side=1, two.sided=TRUE)
add.log.axis(logs=c(1),   tcl=+1.3*abs(par()$tcl), side=1, two.sided=TRUE)
add.log.axis(logs=c(1, 2, 3, 4, 6), side=2, make.labs=TRUE, las=1,
             label="National Hydrography Dataset Basin Area, in km^2")
add.log.axis(logs=c(1, 2, 3, 4, 6), side=1, make.labs=TRUE, las=1,
             label="USGS National Water Information System Contributing Drainage Area, in km^2")
abline(0,1)
abline(1/3,1,lty=2); abline(-1/3,1,lty=2)
abline(1/2,1,lty=2); abline(-1/2,1,lty=3)
mtext("Diagnostic check on watershed areas")
jnk <- abs(DDo$basin_area - DDo$CDA)
summary(jnk[jnk > 1/2])
sites_of_area_bust <- unique(DDo$site_no[jnk > 1/2])
DD_sites_of_area_bust <- DDo[DDo$site_no == sites_of_area_bust[1], ]
for(site in sites_of_area_bust[2:length(sites_of_area_bust)]) {
  DD_sites_of_area_bust <- rbind(DD_sites_of_area_bust,DDo[DDo$site_no == site, ] )
}
DD <- DDo; bust <- 0; k <- 0
for(site in sites_of_area_bust) {
  k <- k + 1
  points(10^DDo$CDA[DDo$site_no == site], 10^DDo$basin_area[DDo$site_no == site], pch=16, col=2)
  bust[k] <- abs(DD$basin_area[DD$site_no == site] - DD$CDA[DD$site_no == site])[1]
  DD <- DD[DD$site_no != site,]
}
text(0,5, paste(sites_of_area_bust, collapse=", "), cex=0.6, pos=4)
par(mgp=c(3,1,0)) # restore defaults
dev.off()


message("Mean bust: ", round(mean(bust), digits=2))
message("Max bust: ",  round(max( bust), digits=2))

print(summary(DDo$CDA - DDo$basin_area))
print(summary( DD$CDA -  DD$basin_area))

#"points(DD$CDA[DD$site_no == "08167000"],
#       DD$basin_area[DD$site_no == "08167000"], pch=16, col=4)

DD$site_no[DD$nid_storage < DD$norm_storage]
#[1] "02295420" "02296750"
DD$decade[DD$nid_storage < DD$norm_storage]
#[1] 2000 2000


ks <- data.frame(decade=DD$decade, cex=NA)
k <- 0; kss <- c(0.6,0.8,1.0,1.2,1.4,1.6)
for(d in sort(unique(ks$decade))) {
   k <- k + 1
   ks$cex[ks$decade == d] <- kss[k]
}
pdf("second_diagnostic_on_area.pdf", useDingbats=FALSE, width=7.5, height=6.5)
par(mgp=c(3,0.5,0)) # going to tick to the inside, change some parameters
plot(10^DD$CDA, 10^DD$basin_area, log="xy", pch=1, lwd=1,
     xaxt="n", yaxt="n", cex=ks$cex, xlim=c(1,300000), ylim=c(1,300000),
     col=rgb(1-as.numeric(DD$pplo == 0),0,as.numeric(DD$pplo == 0),.5),
     xlab="", ylab="")
abline(0,1)
add.log.axis(side=2,      tcl=0.8*abs(par()$tcl), two.sided=TRUE)
add.log.axis(logs=c(1),   tcl=-0.5*abs(par()$tcl), side=2, two.sided=TRUE)
add.log.axis(logs=c(1),   tcl=+1.3*abs(par()$tcl), side=2, two.sided=TRUE)
add.log.axis(side=1,      tcl=0.8*abs(par()$tcl), two.sided=TRUE)
add.log.axis(logs=c(1),   tcl=-0.5*abs(par()$tcl), side=1, two.sided=TRUE)
add.log.axis(logs=c(1),   tcl=+1.3*abs(par()$tcl), side=1, two.sided=TRUE)
add.log.axis(logs=c(1, 2, 3, 4, 6), side=2, make.labs=TRUE, las=1,
             label="National Hydrography Dataset basin area, in km^2")
add.log.axis(logs=c(1, 2, 3, 4, 6), side=1, make.labs=TRUE, las=1,
             label="USGS National Water Information System Contributing Drainage Area, in km^2")

legend(1,200000, c("Watershed for which no zero-flow days were observed (size is decade)",
                   "Watershed for which some zero-flow days were observed (size is decade)"), bty="n",
          col=c(4,2), pch=c(1,1), cex = 0.8)
par(mgp=c(3,1,0)) # restore defaults
dev.off()

D <- DD;

D <- D[D$ed_rch_zone != "1",]

D <- D[D$site_no != "02341500", ]
# Streamgages 02341500 and 02341505 described in the previous section are
# nearly duplicate records and share COMID 3434284 and HUC12 031300030104.
# Only streamgage 02341505, having one extra decade, was retained for
# construction of statistical models for this study.


duan_smearing_estimator <- function(model) { sum(10^residuals(model))/length(residuals(model)) }

save(bnd, D, DD, DDo, knots,
     DD_sites_of_area_bust, duan_smearing_estimator, file="FDCEST.RData")

#  [1] "site_no"             "comid"               "huc12"               "decade"              "dec_long_va"
#  [6] "dec_lat_va"          "n"                   "nzero"               "pplo"                "min"
# [11] "f0.02"               "f0.05"               "f0.1"                "f0.2"                "f0.5"
# [16] "f01"                 "f02"                 "f05"                 "f10"                 "f20"
# [21] "f25"                 "f30"                 "f40"                 "f50"                 "f60"
# [26] "f70"                 "f75"                 "f80"                 "f90"                 "f95"
# [31] "f98"                 "f99"                 "f99.5"               "f99.8"               "f99.9"
# [36] "f99.95"              "f99.98"              "max"                 "L1"                  "L2"
# [41] "T3"                  "T4"                  "T5"                  "T6"                  "T7"
# [46] "T8"                  "median_nonzero"      "major"               "ndams"               "nid_storage"
# [51] "norm_storage"        "ppt_mean"            "ppt_sd"              "temp_mean"           "temp_sd"
# [56] "barren"              "cultivated_cropland" "deciduous_forest"    "developed"           "evergreen_forest"
# [61] "grassland"           "hay_pasture"         "herbaceous_wetland"  "mixed_forest"        "shrubland"
# [66] "water"               "woody_wetland"       "bfi"                 "sinuosity"           "length_km"
# [71] "strm_dens"           "twi"                 "basin_area"          "cat_soller"          "soller"
# [76] "cat_aquifers"        "aquifers"            "bedperm"             "cat_physio"          "physio"
# [81] "cat_ecol3"           "ecol3"               "hlr"                 "rdx"                 "basin_slope"
# [86] "elev_mean"           "statsgo"             "flood_storage"       "ed_rch_zone"         "comid.y"
# [91] "site_no.y"           "huc12.y"             "dec_long_va.y"       "dec_lat_va.y"        "decade.y"
# [96] "dni_ann"             "dni_jan"             "dni_feb"             "dni_mar"             "dni_apr"
#[101] "dni_may"             "dni_jun"             "dni_jul"             "dni_aug"             "dni_sep"
#[106] "dni_oct"             "dni_nov"             "dni_dec"             "CDA"                 "east"
#[111] "north"               "x"                   "y"

# It appears critical that the boundary have variables named say x,y
# The knots have the same names (x,y) and most difficult to figure out
# the x, y must be passed as same names in to the s(...). The x,y is not
# the critical piece, it is that they are all the same literal string.
# For example, v,w would work too.
plot(bnd[[1]]$x,bnd[[1]]$y,type="l", col=8, lwd=.6)
points(D$east[D$nzero == 0], D$north[D$nzero == 0], pch=4, lwd=.5, cex=0.9, col=rgb(0,.5,0.5,.5))
points(D$east[D$nzero > 0 & D$nzero <= 800], D$north[D$nzero > 0 & D$nzero <= 800], pch=4, lwd=.5, cex=0.9, col=rgb(1,0,0.5,.5))
points(D$east[D$nzero > 3653-2000], D$north[D$nzero > 3653-2000], pch=16, lwd=.5, cex=0.9, col=rgb(0.5,0,1))


pdf("test_survreg_gam.pdf", useDingbats=FALSE, height=6, width=7.5)
par(las=1)
family <- "gaussian"
Z <- D
Z$flowtime <- log10(Z$n - Z$nzero)
Zc <- Surv(Z$flowtime, Z$nzero != 0, type="right")
SMt <- survreg(Zc~basin_area+ppt_mean+decade-1, data=Z, dist=family)
Pt <- Pto <- predict(SMt); Pt[Pt > log10(3653)] <- log10(3653)

Z <- D
x <- Z$east; y <- Z$north
Z$left.threshold <-  log10(rep(0, length(Z$nzero)))
Z$right.threshold <- log10(Z$n)
Z$flowtime <- log10(Z$n - Z$nzero)
GMt <- gam(flowtime~basin_area+s(ppt_mean, k=2)+decade-1,
           family=tobit1(left.threshold=  Z$left.threshold,
                         right.threshold=Z$right.threshold), data=Z); summary(GMt)
Gt <- Gto <- predict(GMt); Gt[Gt > log10(3653)] <- log10(3653)
plot(Gto, Pto, ylab="Test survival regression: log10(days of decadal streamflow)",
               xlab="Test censored generalized additive model (GAM): log10(days of decadal streamflow)",
               pch=21, col=4, bg=8, cex=0.5, lwd=0.5, xaxs="i", yaxs="i",
               xlim=c(3.3,4.0), ylim=c(3.3,4.0), type="n")
usr <- par()$usr
xp <- c(log10(3653),  4, 4, 3.3, 3.3, log10(3653), log10(3653))
yp <- c(3.3, 3.3, 4, 4, log10(3653), log10(3653), 3.3)

polygon(xp, yp, col=grey(0.97), lty=0)
ix <- sample(1:length(Gto)) # randomizer
points(Gto[ix], Pto[ix], pch=21, col=grey(0.5), bg=Z$decade[ix], cex=0.7, lwd=0.5)
lines(rep(log10(3653),2), c(3.3, 4), lty=2, lwd=0.65, col=grey(0.5))
lines(c(3.3, 4), rep(log10(3653),2), lty=2, lwd=0.65, col=grey(0.5))

abline(0,1, lwd=0.8)
length(Z$n[Z$nzero == 0])
#[1] 2002
length(Z$n[Z$nzero != 0])
#[1] 739
length(Pto[Pto < log10(3653)])
#[1] 387
length(Gto[Gto < log10(3653)])
#[1] 414
txt <- paste0("Number of streamgage:decade records with no zero-flow conditions: 2,002\n",
              "Number of streamgage:decade records with zero-flow conditions: 739\n",
              "Number of survival regression predicted with zero-flow conditions: 387\n",
              "Number of GAM regression predicted with zero-flow conditions: 414\n")
text(3.57, 3.4, txt, pos=4, cex=0.7)
txt <- paste0("All three grey regions depict predictions of\n",
              "no zero-flow occurrences in at least a decade.\n",
              "These exemplify the concept of right censoring.\n",
              "Given a time period longer than a decade, one\n",
              "or more zero-flow might have been observed.")
text(3.295, 3.9, txt, pos=4, cex=0.7)
text(3.95, 3.95, "Equal value line", cex=0.7)
legend(3.31, 3.75, c("A unique streamgage:decade record\n in the database colored by decade"),
              pch=21, lty=0, col=grey(0.5), pt.bg=2, bty="n", cex=0.8)
dev.off()



family <- "gaussian"
Z <- D
x <- Z$east; y <- Z$north
Z$flowtime <- log10(Z$n - Z$nzero)
Zc <- Surv(Z$flowtime, Z$nzero != 0, type="right")
SM0 <- survreg(Zc~basin_area+
                    ppt_mean+temp_mean+dni_ann+
                    developed+grassland+
                    bedperm+decade-1, data=Z, dist=family)
Psurv0 <- Psurv0o <- predict(SM0); Psurv0[Psurv0 > log10(3653)] <- log10(3653)

Z <- D
x <- Z$east; y <- Z$north
Z$left.threshold <-  log10(rep(0, length(Z$nzero)))
Z$right.threshold <- log10(Z$n)
Z$flowtime <- log10(Z$n - Z$nzero)
GM0 <- gam(flowtime~basin_area+
                    ppt_mean+s(temp_mean)+dni_ann+
                    developed+s(grassland)+
                    bedperm+decade-1,
           family=tobit1(left.threshold=  Z$left.threshold,
                         right.threshold=Z$right.threshold), data=Z); summary(GM0)
P0 <- P0o <- predict(GM0); P0[P0 > log10(3653)] <- log10(3653)
plot(Psurv0, P0, xlab="Survival Regression: log10(days of decadal streamflow)",
                 ylab="Censored GAM: log10(days of decadal streamflow)")
abline(0,1)
# August 15, 2018 testing with the "final dataset" does not produce convergence warnings.
#mtext("Model structure not quite exact because of convergence issues in GAM")
# WHA I would rather not have the smooth on temp_mean and grassland but these
# appear needed to get the GAM to converge. So the model does not quite match
# that from survreg().

GM1 <- gam(flowtime~s(basin_area)+
                    s(ppt_mean, k=5)+s(temp_mean, k=4)+s(dni_ann, k=7)+
                    developed+grassland+
                    bedperm+decade-1,
           family=tobit1(left.threshold=  Z$left.threshold,
                         right.threshold=Z$right.threshold), data=Z); summary(GM1)
GM2 <- gam(flowtime~s(basin_area)+
                    s(ppt_mean, k=5)+s(temp_mean, k=4)+s(dni_ann, k=7)+
                    developed+grassland+
                    bedperm+decade-1+
        s(x,y, bs="so", xt=list(bnd=bnd)), knots=knots,
           family=tobit1(left.threshold=  Z$left.threshold,
                         right.threshold=Z$right.threshold), data=Z); summary(GM2)
GM3 <- gam(flowtime~s(basin_area, k=5)+
                    s(ppt_mean, k=5)+s(temp_mean, k=4)+s(dni_ann, k=7)+
                    developed+grassland+
                    bedperm+decade-1+
        s(x,y), knots=knots,
           family=tobit1(left.threshold=  Z$left.threshold,
                         right.threshold=Z$right.threshold), data=Z); summary(GM3)
GM4 <- gam(flowtime~s(basin_area, k=5)+
                    s(ppt_mean, k=5)+s(temp_mean, k=4)+s(dni_ann, k=7)+s(dni_mar, k=7)+
                    developed+grassland+
                    bedperm+decade-1+
        s(x,y), knots=knots,
           family=tobit1(left.threshold=  Z$left.threshold,
                         right.threshold=Z$right.threshold), data=Z); summary(GM4)

P <- Po <- predict(SM0); P[P > log10(3653)] <- log10(3653)
C1 <- C1o <- predict(GM1); C1[C1 > log10(3653)] <- log10(3653)
C2 <- C2o <- predict(GM2); C2[C2 > log10(3653)] <- log10(3653)
C3 <- C3o <- predict(GM3); C3[C3 > log10(3653)] <- log10(3653)
C4 <- C4o <- predict(GM4); C4[C4 > log10(3653)] <- log10(3653)
pa <- abs(P - Z$flowtime)
ca <- abs(C1- Z$flowtime)
cb <- abs(C2- Z$flowtime)
cc <- abs(C3- Z$flowtime)
cd <- abs(C4- Z$flowtime)
summary(pa)
summary(ca)
summary(cb)
summary(cc)
summary(cd)


plot(Z$flowtime, P, xlim=c(2.9,3.6), ylim=c(3,3.6), type="n")
abline(0,1)
points(Z$flowtime ,P, pch=16, col=rgb(0,0,1,.5), cex=0.5)
points(Z$flowtime,C2, pch=16, col=rgb(1,0,0,.2), cex=0.5)
for(i in 1:length(P)) {
  if(abs(P[i]-C2[i]) < .02) next
  try(arrows(Z$flowtime[i],P[i],Z$flowtime[i],C2[i],
             lwd=0.5, angle=20, length=.1, col=2), silent=TRUE)
}

plot(  qnorm(pp(pa)), sort(pa), type="l", xlim=c(0,4))
lines( qnorm(pp(ca)), sort(ca), col=4)
lines( qnorm(pp(cb)), sort(cb), col=2)
lines( qnorm(pp(cc)), sort(cc), col=3)
lines( qnorm(pp(cd)), sort(cd), col=6)
points(qnorm(pp(cd)), sort(cd), col=6, lwd=0.5)

Ppplo  <- (3653-10^P) /3653
C1pplo <- (3653-10^C1)/3653
C2pplo <- (3653-10^C2)/3653
C3pplo <- (3653-10^C3)/3653
plot(  Z$pplo, Ppplo, xlim=c(0,1), ylim=c(0,1))
points(Z$pplo, C2pplo, col=2)
points(Z$pplo, C2pplo, col=4)
points(Z$pplo, C3pplo, col=3)
abline(0,1)

PPLO <- GM4
save(DDo, DD, D, SM0, GM1, GM2, GM3, GM4, PPLO, file="PPLOS.RData")

pplo.sigma <- PPLO$pplo.sigma <- sqrt(mean((predict(PPLO) - Z$flowtime)^2))
PGAM <- gamIntervals(predict(PPLO, se.fit=TRUE), gam=PPLO, interval="prediction", sigma=pplo.sigma)
# In an uncensored data world, the following will be a 1:1 relation (see gamIntervals), but we won't fully
# see that when there is the censoring. But the abline will plot through the middle of the data cloud.
#plot(PPLO$hat, (PGAM$se.fit/PGAM$residual.scale)^2)
#abline(0,1)

# Terms invert in upper/lower meaning and hence the flipping during data.frame construction.
PPLOdf <- data.frame(comid=Z$comid, site_no=Z$site_no, huc12=Z$huc12,
                     decade=Z$decade,
                     dec_long_va=Z$dec_long_va, dec_lat_va=Z$dec_lat_va,
                     in_model_pplo=TRUE, pplo=Z$pplo, flowtime=log10(Z$n - Z$nzero),
                     est_lwr_pplo=(3653-10^PGAM$upr)/3653,
                     est_pplo    =(3653-10^PGAM$fit)/3653,
                     est_upr_pplo=(3653-10^PGAM$lwr)/3653,
                     est_lwr_flowtime=PGAM$lwr,
                     est_flowtime=PGAM$fit,
                     est_upr_flowtime=PGAM$upr,
                     stringsAsFactors=FALSE)
PPLOdf$est_lwr_pplo[PPLOdf$est_lwr_pplo < 0] <- 0
PPLOdf$est_pplo[    PPLOdf$est_pplo     < 0] <- 0
PPLOdf$est_upr_pplo[PPLOdf$est_upr_pplo < 0] <- 0
PPLOdf$rse_pplo <- pplo.sigma
PPLOdf$se.fit_pplo <- PGAM$se.fit

sites_to_fill <- unique(c(sites_of_area_bust, DDo$site_no[DDo$ed_rch_zone == 1]))
i <- 0
for(site in sites_to_fill) {
  for(decade in as.character(DDo$decade[DDo$site_no == site])) {
    i <- i + 1
    message(site," ", decade, " ", i)
    tmp <- data.frame(comid=DDo$comid[DDo$site_no == site & DDo$decade == decade],
                      site_no=site,
                      huc12=DDo$huc12[DDo$site_no == site & DDo$decade == decade],
                      decade=decade,
                      dec_long_va=DDo$dec_long_va[DDo$site_no == site & DDo$decade == decade],
                      dec_lat_va=DDo$dec_lat_va[DDo$site_no == site & DDo$decade == decade],
                      in_model_pplo=FALSE,
                      pplo=DDo$pplo[DDo$site_no == site & DDo$decade == decade],
                      flowtime=log10(DDo$n[    DDo$site_no == site & DDo$decade == decade] -
                                     DDo$nzero[DDo$site_no == site & DDo$decade == decade]),
                      est_lwr_pplo=NA, est_pplo=NA, est_upr_pplo=NA,
                      est_lwr_flowtime=NA, est_flowtime=NA, est_upr_flowtime=NA,
                      rse_pplo=NA, se.fit_pplo=NA, stringsAsFactors=FALSE)
    PPLOdf <- rbind(PPLOdf, tmp)
  }
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(PPLO, newdata=tmp, se.fit=TRUE)
}
PPLOdf <- PPLOdf[order(PPLOdf$site_no, PPLOdf$decade),]

for(site in sites_to_fill) {
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(PPLO, newdata=tmp, se.fit=TRUE)
  pgk <- gamIntervals(jnk, gam=PPLO, interval="prediction", sigma=pplo.sigma)
  df <- data.frame(comid=tmp$comid, site_no=tmp$site_no, huc12=tmp$huc12,
                   decade=tmp$decade,
                   dec_long_va=tmp$dec_long_va, dec_lat_va=tmp$dec_lat_va,
                   in_model_pplo=0,
                   pplo=tmp$pplo, flowtime=log10(tmp$n - tmp$nzero),
                   est_lwr_pplo=(3653-10^pgk$upr)/3653,
                   est_pplo    =(3653-10^pgk$fit)/3653,
                   est_upr_pplo=(3653-10^pgk$lwr)/3653,
                   est_lwr_flowtime=pgk$lwr,
                   est_flowtime=pgk$fit,
                   est_upr_flowtime=pgk$upr,
                   stringsAsFactors=FALSE)
  df$est_lwr_pplo[df$est_lwr_pplo < 0] <- 0
  df$est_pplo[    df$est_pplo     < 0] <- 0
  df$est_upr_pplo[df$est_upr_pplo < 0] <- 0
  df$rse_pplo <- pplo.sigma
  df$se.fit_pplo <- pgk$se.fit
  PPLOdf[PPLOdf$site_no == site,] <- df
}
#PPLOdf$in_model_pplo[PPLOdf$in_model_pplo == 1] <- "yes"
#PPLOdf$in_model_pplo[PPLOdf$in_model_pplo == 0] <- "no"
summary(PPLOdf$pplo <= PPLOdf$est_upr_pplo*0.62) # in bulk the limits just are wrong for the
# pplo but that is looking at them as if they were noncensored that is compounded by the
# conversion of flowtime to a pplo (probability). If we have 2804 rows, and if we continue
# to expect the upper limit to be a 0.025 probability exceeding, then we expect 70 cases
# as predictions outside. So this is a rough fix. Will need to evaluate this in cross-validation.
#$PPLOdf$est_upr_pplo <- 0.62*PPLOdf$est_upr_pplo
# All of the est_lwr_pplo are zero anyway so know simple way seen how to fix them.

sum(abs(PPLOdf$est_pplo - PPLOdf$pplo) <= 0   )/length(DDo$site_no)
sum(abs(PPLOdf$est_pplo - PPLOdf$pplo) <= 0.02)/length(DDo$site_no)
sum(abs(PPLOdf$est_pplo - PPLOdf$pplo) <= 0.05)/length(DDo$site_no)
sum(abs(PPLOdf$est_pplo - PPLOdf$pplo) <= 0.10)/length(DDo$site_no)

#write_feather(PPLOdf, "all_gage_est_pplo.feather")



plot(bnd[[1]]$x,bnd[[1]]$y,type="l", col=8, lwd=.6)
points(D$east[D$nzero == 0], D$north[D$nzero == 0], pch=4, lwd=.5, cex=0.9, col=rgb(0,.5,0.5,.5))
points(D$east[D$nzero > 0 & D$nzero <= 800], D$north[D$nzero > 0 & D$nzero <= 800], pch=4, lwd=.5, cex=0.9, col=rgb(1,0,0.5,.5))
points(D$east[D$nzero > 3653-2000], D$north[D$nzero > 3653-2000], pch=16, lwd=.5, cex=0.9, col=rgb(0.5,0,1))
mtext("Raw of PPLO")

plot(bnd[[1]]$x,bnd[[1]]$y,type="l", col=8, lwd=.6)
points(D$east[P == log10(3653)], D$north[P == log10(3653)], pch=4, lwd=.5, cex=0.9, col=rgb(0,.5,0.5,.5))
points(D$east[P < log10(3653)], D$north[P < log10(3653)], pch=4, lwd=.5, cex=0.9, col=rgb(1,0,0.5,.5))
points(D$east[P < log10(2000)], D$north[P < log10(2000)], pch=16, lwd=.5, cex=0.9, col=rgb(0.5,0,1,.5))
mtext("Survival Regression Predictions of PPLO")

plot(bnd[[1]]$x,bnd[[1]]$y,type="l", col=8, lwd=.6)
points(D$east[C1 == log10(3653)], D$north[C1 == log10(3653)], pch=4, lwd=.5, cex=0.9, col=rgb(0,.5,0.5,.5))
points(D$east[C1 < log10(3653)], D$north[C1 < log10(3653)], pch=4, lwd=.5, cex=0.9, col=rgb(1,0,0.5,.5))
points(D$east[C1 < log10(2000)], D$north[C1 < log10(2000)], pch=16, lwd=.5, cex=0.9, col=rgb(0.5,0,1,.5))
mtext("Predictions of PPLO (C1)")

plot(bnd[[1]]$x,bnd[[1]]$y,type="l", col=8, lwd=.6)
points(D$east[C2 == log10(3653)], D$north[C2 == log10(3653)], pch=4, lwd=.5, cex=0.9, col=rgb(0,.5,0.5,.5))
points(D$east[C2 < log10(3653)], D$north[C2 < log10(3653)], pch=4, lwd=.5, cex=0.9, col=rgb(1,0,0.5,.5))
points(D$east[C2 < log10(2000)], D$north[C2 < log10(2000)], pch=16, lwd=.5, cex=0.9, col=rgb(0.5,0,1,.5))
mtext("Predictions cenGAM of PPLO (C2)")

plot(bnd[[1]]$x,bnd[[1]]$y,type="l", col=8, lwd=.6)
points(D$east[C3 == log10(3653)], D$north[C3 == log10(3653)], pch=4, lwd=.5, cex=0.9, col=rgb(0,.5,0.5,.5))
points(D$east[C3 < log10(3653)], D$north[C3 < log10(3653)], pch=4, lwd=.5, cex=0.9, col=rgb(1,0,0.5,.5))
points(D$east[C3 < log10(2000)], D$north[C3 < log10(2000)], pch=16, lwd=.5, cex=0.9, col=rgb(0.5,0,1,.5))
mtext("Predictions cenGAM of PPLO (C3)")

pdf("PPLO.pdf", useDingbats=FALSE)
  plot(Z$flowtime, fitted.values(PPLO), col=3)
  abline(0,1)
  plot(PPLO, scheme=2, residuals=TRUE)
  points(Z$x, Z$y, pch=4, lwd=.5, cex=0.9, col=8)
  points(knots$x, knots$y, pch=16, cex=1.1, col=4)
  text(100, 500, "PPLO")
dev.off()

###### END PPLO


Z <- D
x <- Z$x; y <- Z$y
z <- log10(Z$L1)
L1   <- gam(z~s(basin_area) + s(basin_slope, k=5)+
              s(ppt_mean, k=5)+ s(dni_ann, k=7)+
              developed+
              decade-1+
              s(x, y), knots=knots, data=Z, # , bs="so", xt=list(bnd=bnd)
              family="gaussian")
L1$duan_smearing <- duan_smearing_estimator(L1)
pdf("L1.pdf", useDingbats=FALSE)
  plot(z, fitted.values(L1), col=3)
  abline(0,1)
  plot(L1, scheme=2, residuals=TRUE)
  points(Z$x, Z$y, pch=4, lwd=.5, cex=0.9, col=8)
  points(knots$x, knots$y, pch=16, cex=1.1, col=4)
  text(100, 500, "L1")
dev.off()


PGAM <- gamIntervals(predict(L1, se.fit=TRUE), gam=L1, interval="prediction")
#plot(L1$hat, (PGAM$se.fit/PGAM$residual.scale)^2)
#abline(0,1)
sigma <- PGAM$residual.scale[1]
L1df <- data.frame(comid=Z$comid, site_no=Z$site_no, huc12=Z$huc12,
                     decade=Z$decade,
                     dec_long_va=Z$dec_long_va, dec_lat_va=Z$dec_lat_va,
                     in_model_L1=TRUE, bias_corr=L1$duan_smearing, L1=Z$L1,
                     est_lwr_L1=10^PGAM$lwr,
                     est_L1    =10^PGAM$fit,
                     est_upr_L1=10^PGAM$upr, stringsAsFactors=FALSE)
L1df$rse_L1 <- sigma
L1df$se.fit_L1 <- PGAM$se.fit

sites_to_fill <- unique(c(sites_of_area_bust, DDo$site_no[DDo$ed_rch_zone == 1]))
i <- 0
for(site in sites_to_fill) {
  for(decade in as.character(DDo$decade[DDo$site_no == site])) {
    i <- i + 1
    message(site," ", decade, " ", i)
    tmp <- data.frame(comid=DDo$comid[DDo$site_no == site & DDo$decade == decade],
                      site_no=site,
                      huc12=DDo$huc12[DDo$site_no == site & DDo$decade == decade],
                      decade=decade,
                      dec_long_va=DDo$dec_long_va[DDo$site_no == site & DDo$decade == decade],
                      dec_lat_va=DDo$dec_lat_va[DDo$site_no == site & DDo$decade == decade],
                      in_model_L1=FALSE, bias_corr=L1$duan_smearing,
                      L1=DDo$L1[DDo$site_no == site & DDo$decade == decade],
                      est_lwr_L1=NA, est_L1=NA, est_upr_L1=NA,
                      rse_L1=NA, se.fit_L1=NA, stringsAsFactors=FALSE)
    L1df <- rbind(L1df, tmp)
  }
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(L1, newdata=tmp, se.fit=TRUE)
}
L1df <- L1df[order(L1df$site_no, L1df$decade),]

for(site in sites_to_fill) {
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(L1, newdata=tmp, se.fit=TRUE)
  pgk <- gamIntervals(jnk, gam=L1, interval="prediction")
  df <- data.frame(comid=tmp$comid, site_no=tmp$site_no, huc12=tmp$huc12,
                   decade=tmp$decade,
                   dec_long_va=tmp$dec_long_va, dec_lat_va=tmp$dec_lat_va,
                   in_model_L1=0, bias_corr=L1$duan_smearing,
                   L1=tmp$L1,
                   est_lwr_L1=10^pgk$lwr,
                   est_L1    =10^pgk$fit,
                   est_upr_L1=10^pgk$upr, stringsAsFactors=FALSE)
  df$rse_L1 <- sigma
  df$se.fit_L1 <- pgk$se.fit
  L1df[L1df$site_no == site,] <- df
}

#write_feather(L1df, "all_gage_est_L1.feather")

z <- D$L2/D$L1 # --------------------------- Coefficient of L-variation
T2   <- gam(z~s(basin_area) + s(basin_slope, k=5)+
              s(ppt_mean, k=5) + s(temp_mean, k=4) + s(dni_ann, k=7)+
              developed+
              s(flood_storage, k=5)+
              decade-1+
              s(x, y), knots=knots, data=D, # , bs="so", xt=list(bnd=bnd)
              family="gaussian")
pdf("T2.pdf", useDingbats=FALSE)
  plot(z, fitted.values(T2))
  abline(0,1)
  plot(T2, scheme=2, residuals=TRUE)

  points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
  points(knots$x, knots$y, pch=16, cex=1.1, col=4)
  text(100, 500, "Tau2")
dev.off()



PGAM <- gamIntervals(predict(T2, se.fit=TRUE), gam=T2, interval="prediction")
#plot(T2$hat, (PGAM$se.fit/PGAM$residual.scale)^2)
#abline(0,1)
sigma <- PGAM$residual.scale[1]
# Terms invert in upper/lower meaning and hence the flipping during data.frame construction.
T2df <- data.frame(comid=Z$comid, site_no=Z$site_no, huc12=Z$huc12,
                     decade=Z$decade,
                     dec_long_va=Z$dec_long_va, dec_lat_va=Z$dec_lat_va,
                     in_model_T2=TRUE, T2=Z$L2/Z$L1,
                     est_lwr_T2=PGAM$lwr,
                     est_T2    =PGAM$fit,
                     est_upr_T2=PGAM$upr, stringsAsFactors=FALSE)
T2df$rse_T2 <- sigma
T2df$se.fit_T2 <- PGAM$se.fit

sites_to_fill <- unique(c(sites_of_area_bust, DDo$site_no[DDo$ed_rch_zone == 1]))
i <- 0
for(site in sites_to_fill) {
  for(decade in as.character(DDo$decade[DDo$site_no == site])) {
    i <- i + 1
    message(site," ", decade, " ", i)
    tmp <- data.frame(comid=DDo$comid[DDo$site_no == site & DDo$decade == decade],
                      site_no=site,
                      huc12=DDo$huc12[DDo$site_no == site & DDo$decade == decade],
                      decade=decade,
                      dec_long_va=DDo$dec_long_va[DDo$site_no == site & DDo$decade == decade],
                      dec_lat_va=DDo$dec_lat_va[DDo$site_no == site & DDo$decade == decade],
                      in_model_T2=FALSE,
                      T2=DDo$L2[DDo$site_no == site & DDo$decade == decade]/DDo$L1[DDo$site_no == site & DDo$decade == decade],
                      est_lwr_T2=NA, est_T2=NA, est_upr_T2=NA,
                      rse_T2=NA, se.fit_T2=NA, stringsAsFactors=FALSE)
    T2df <- rbind(T2df, tmp)
  }
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(T2, newdata=tmp, se.fit=TRUE)
}
T2df <- T2df[order(T2df$site_no, T2df$decade),]

for(site in sites_to_fill) {
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(T2, newdata=tmp, se.fit=TRUE)
  pgk <- gamIntervals(jnk, gam=T2, interval="prediction")
  df <- data.frame(comid=tmp$comid, site_no=tmp$site_no, huc12=tmp$huc12,
                   decade=tmp$decade,
                   dec_long_va=tmp$dec_long_va, dec_lat_va=tmp$dec_lat_va,
                   in_model_T2=0,
                   T2=tmp$L2/tmp$L1,
                   est_lwr_T2=pgk$lwr,
                   est_T2    =pgk$fit,
                   est_upr_T2=pgk$upr, stringsAsFactors=FALSE)
  df$rse_T2 <- sigma
  df$se.fit_T2 <- pgk$se.fit
  T2df[T2df$site_no == site,] <- df
}

#write_feather(T2df, "all_gage_est_T2.feather")




z <- D$T3      # --------------------------- L-skew
T3   <- gam(z~s(basin_area)+ s(basin_slope, k=5)+
              s(ppt_mean, k=5) + s(temp_mean, k=4) + s(dni_ann, k=7)+
              s(flood_storage, k=5)+
              decade-1+
              s(x, y), knots=knots, data=D, #, bs="so", xt=list(bnd=bnd)
              family="gaussian")
pdf("T3.pdf", useDingbats=FALSE)
  plot(z, fitted.values(T3))
  abline(0,1)
  plot(T3, scheme=2, residuals=TRUE)
  points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
  points(knots$x, knots$y, pch=16, cex=1.1, col=4)
  text(100, 500, "Tau3")
dev.off()


PGAM <- gamIntervals(predict(T3, se.fit=TRUE), gam=T3, interval="prediction")
#plot(T3$hat, (PGAM$se.fit/PGAM$residual.scale)^2)
#abline(0,1)
sigma <- PGAM$residual.scale[1]
# Terms invert in upper/lower meaning and hence the flipping during data.frame construction.
T3df <- data.frame(comid=Z$comid, site_no=Z$site_no, huc12=Z$huc12,
                     decade=Z$decade,
                     dec_long_va=Z$dec_long_va, dec_lat_va=Z$dec_lat_va,
                     in_model_T3=TRUE, T3=Z$T3,
                     est_lwr_T3=PGAM$lwr,
                     est_T3    =PGAM$fit,
                     est_upr_T3=PGAM$upr, stringsAsFactors=FALSE)
T3df$rse_T3 <- sigma
T3df$se.fit_T3 <- PGAM$se.fit

sites_to_fill <- unique(c(sites_of_area_bust, DDo$site_no[DDo$ed_rch_zone == 1]))
i <- 0
for(site in sites_to_fill) {
  for(decade in as.character(DDo$decade[DDo$site_no == site])) {
    i <- i + 1
    message(site," ", decade, " ", i)
    tmp <- data.frame(comid=DDo$comid[DDo$site_no == site & DDo$decade == decade],
                      site_no=site,
                      huc12=DDo$huc12[DDo$site_no == site & DDo$decade == decade],
                      decade=decade,
                      dec_long_va=DDo$dec_long_va[DDo$site_no == site & DDo$decade == decade],
                      dec_lat_va=DDo$dec_lat_va[DDo$site_no == site & DDo$decade == decade],
                      in_model_T3=FALSE,
                      T3=DDo$T3[DDo$site_no == site & DDo$decade == decade],
                      est_lwr_T3=NA, est_T3=NA, est_upr_T3=NA,
                      rse_T3=NA, se.fit_T3=NA, stringsAsFactors=FALSE)
    T3df <- rbind(T3df, tmp)
  }
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(T3, newdata=tmp, se.fit=TRUE)
}
T3df <- T3df[order(T3df$site_no, T3df$decade),]

for(site in sites_to_fill) {
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(T3, newdata=tmp, se.fit=TRUE)
  pgk <- gamIntervals(jnk, gam=T3, interval="prediction")
  df <- data.frame(comid=tmp$comid, site_no=tmp$site_no, huc12=tmp$huc12,
                   decade=tmp$decade,
                   dec_long_va=tmp$dec_long_va, dec_lat_va=tmp$dec_lat_va,
                   in_model_T3=0,
                   T3=tmp$T3,
                   est_lwr_T3=pgk$lwr,
                   est_T3    =pgk$fit,
                   est_upr_T3=pgk$upr, stringsAsFactors=FALSE)
  df$rse_T3 <- sigma
  df$se.fit_T3 <- pgk$se.fit
  T3df[T3df$site_no == site,] <- df
}

#write_feather(T3df, "all_gage_est_T3.feather")





z <- D$T4      # --------------------------- L-kurtosis
T4   <- gam(z~s(basin_area) +s(basin_slope, k=5)+
              s(ppt_mean, k=5) + s(temp_mean, k=4) + s(dni_ann, k=7)+
              s(flood_storage, k=5)+
              decade-1+
              s(x, y), knots=knots, data=D, # , bs="so", xt=list(bnd=bnd)
              family="gaussian")
pdf("T4.pdf", useDingbats=FALSE)
  plot(z, fitted.values(T4))
  abline(0,1)
  plot(T4, scheme=2, residuals=TRUE)
  points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
  points(knots$x, knots$y, pch=16, cex=1.1, col=4)
  text(100, 500, "Tau4")
dev.off()


PGAM <- gamIntervals(predict(T4, se.fit=TRUE), gam=T4, interval="prediction")
#plot(T4$hat, (PGAM$se.fit/PGAM$residual.scale)^2)
#abline(0,1)
sigma <- PGAM$residual.scale[1]
# Terms invert in upper/lower meaning and hence the flipping during data.frame construction.
T4df <- data.frame(comid=Z$comid, site_no=Z$site_no, huc12=Z$huc12,
                     decade=Z$decade,
                     dec_long_va=Z$dec_long_va, dec_lat_va=Z$dec_lat_va,
                     in_model_T4=TRUE, T4=Z$T4,
                     est_lwr_T4=PGAM$lwr,
                     est_T4    =PGAM$fit,
                     est_upr_T4=PGAM$upr, stringsAsFactors=FALSE)
T4df$rse_T4 <- sigma
T4df$se.fit_T4 <- PGAM$se.fit

sites_to_fill <- unique(c(sites_of_area_bust, DDo$site_no[DDo$ed_rch_zone == 1]))
i <- 0
for(site in sites_to_fill) {
  for(decade in as.character(DDo$decade[DDo$site_no == site])) {
    i <- i + 1
    message(site," ", decade, " ", i)
    tmp <- data.frame(comid=DDo$comid[DDo$site_no == site & DDo$decade == decade],
                      site_no=site,
                      huc12=DDo$huc12[DDo$site_no == site & DDo$decade == decade],
                      decade=decade,
                      dec_long_va=DDo$dec_long_va[DDo$site_no == site & DDo$decade == decade],
                      dec_lat_va=DDo$dec_lat_va[DDo$site_no == site & DDo$decade == decade],
                      in_model_T4=FALSE,
                      T4=DDo$T4[DDo$site_no == site & DDo$decade == decade],
                      est_lwr_T4=NA, est_T4=NA, est_upr_T4=NA,
                      rse_T4=NA, se.fit_T4=NA, stringsAsFactors=FALSE)
    T4df <- rbind(T4df, tmp)
  }
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(T4, newdata=tmp, se.fit=TRUE)
}
T4df <- T4df[order(T4df$site_no, T4df$decade),]

for(site in sites_to_fill) {
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(T4, newdata=tmp, se.fit=TRUE)
  pgk <- gamIntervals(jnk, gam=T4, interval="prediction")
  df <- data.frame(comid=tmp$comid, site_no=tmp$site_no, huc12=tmp$huc12,
                   decade=tmp$decade,
                   dec_long_va=tmp$dec_long_va, dec_lat_va=tmp$dec_lat_va,
                   in_model_T4=0,
                   T4=tmp$T4,
                   est_lwr_T4=pgk$lwr,
                   est_T4    =pgk$fit,
                   est_upr_T4=pgk$upr, stringsAsFactors=FALSE)
  df$rse_T4 <- sigma
  df$se.fit_T4 <- pgk$se.fit
  T4df[T4df$site_no == site,] <- df
}

#write_feather(T4df, "all_gage_est_T4.feather")




z <- D$T5      # --------------------------- Fifth L-moment ratio
T5   <- gam(z~s(basin_area) +s(basin_slope, k=5)+
              s(temp_mean, k=4) + s(dni_ann, k=7)+
              s(flood_storage, k=5)+
              decade-1+
              s(x, y), knots=knots, data=D, #, bs="so", xt=list(bnd=bnd)
              family="gaussian")
pdf("T5.pdf", useDingbats=FALSE)
  plot(z, fitted.values(T5))
  abline(0,1)
  plot(T5, scheme=2, residuals=TRUE)
  points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
  points(knots$x, knots$y, pch=16, cex=1.1, col=4)
  text(100, 500, "Tau5")
dev.off()


PGAM <- gamIntervals(predict(T5, se.fit=TRUE), gam=T5, interval="prediction")
#plot(T5$hat, (PGAM$se.fit/PGAM$residual.scale)^2)
#abline(0,1)
sigma <- PGAM$residual.scale[1]
# Terms invert in upper/lower meaning and hence the flipping during data.frame construction.
T5df <- data.frame(comid=Z$comid, site_no=Z$site_no, huc12=Z$huc12,
                     decade=Z$decade,
                     dec_long_va=Z$dec_long_va, dec_lat_va=Z$dec_lat_va,
                     in_model_T5=TRUE, T5=Z$T5,
                     est_lwr_T5=PGAM$lwr,
                     est_T5    =PGAM$fit,
                     est_upr_T5=PGAM$upr, stringsAsFactors=FALSE)
T5df$rse_T5 <- sigma
T5df$se.fit_T5 <- PGAM$se.fit

sites_to_fill <- unique(c(sites_of_area_bust, DDo$site_no[DDo$ed_rch_zone == 1]))
i <- 0
for(site in sites_to_fill) {
  for(decade in as.character(DDo$decade[DDo$site_no == site])) {
    i <- i + 1
    message(site," ", decade, " ", i)
    tmp <- data.frame(comid=DDo$comid[DDo$site_no == site & DDo$decade == decade],
                      site_no=site,
                      huc12=DDo$huc12[DDo$site_no == site & DDo$decade == decade],
                      decade=decade,
                      dec_long_va=DDo$dec_long_va[DDo$site_no == site & DDo$decade == decade],
                      dec_lat_va=DDo$dec_lat_va[DDo$site_no == site & DDo$decade == decade],
                      in_model_T5=FALSE,
                      T5=DDo$T5[DDo$site_no == site & DDo$decade == decade],
                      est_lwr_T5=NA, est_T5=NA, est_upr_T5=NA,
                      rse_T5=NA, se.fit_T5=NA, stringsAsFactors=FALSE)
    T5df <- rbind(T5df, tmp)
  }
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(T5, newdata=tmp, se.fit=TRUE)
}
T5df <- T5df[order(T5df$site_no, T5df$decade),]

for(site in sites_to_fill) {
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(T5, newdata=tmp, se.fit=TRUE)
  pgk <- gamIntervals(jnk, gam=T5, interval="prediction")
  df <- data.frame(comid=tmp$comid, site_no=tmp$site_no, huc12=tmp$huc12,
                   decade=tmp$decade,
                   dec_long_va=tmp$dec_long_va, dec_lat_va=tmp$dec_lat_va,
                   in_model_T5=0,
                   T5=tmp$T5,
                   est_lwr_T5=pgk$lwr,
                   est_T5    =pgk$fit,
                   est_upr_T5=pgk$upr, stringsAsFactors=FALSE)
  df$rse_T5 <- sigma
  df$se.fit_T5 <- pgk$se.fit
  T5df[T5df$site_no == site,] <- df
}

#write_feather(T5df, "all_gage_est_T5.feather")





z <- D$T6      # --------------------------- Sixth L-moment ratio
T6   <- gam(z~s(basin_area) +
              s(ppt_mean, k=5) + s(temp_mean, k=4) + s(dni_ann, k=7)+
              developed+
              mixed_forest+shrubland+s(flood_storage)+
              decade-1+
              s(x, y), knots=knots, data=D, #, bs="so", xt=list(bnd=bnd)
              family="gaussian")
pdf("T6.pdf", useDingbats=FALSE)
  plot(z, fitted.values(T6))
  abline(0,1)
  plot(T6, scheme=2, residuals=TRUE)
  points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
  points(knots$x, knots$y, pch=16, cex=1.1, col=4)
  text(100, 500, "Tau6")
dev.off()

PGAM <- gamIntervals(predict(T6, se.fit=TRUE), gam=T6, interval="prediction")
#plot(T6$hat, (PGAM$se.fit/PGAM$residual.scale)^2)
#abline(0,1)
sigma <- PGAM$residual.scale[1]
# Terms invert in upper/lower meaning and hence the flipping during data.frame construction.
T6df <- data.frame(comid=Z$comid, site_no=Z$site_no, huc12=Z$huc12,
                     decade=Z$decade,
                     dec_long_va=Z$dec_long_va, dec_lat_va=Z$dec_lat_va,
                     in_model_T6=TRUE, T6=Z$T6,
                     est_lwr_T6=PGAM$lwr,
                     est_T6    =PGAM$fit,
                     est_upr_T6=PGAM$upr, stringsAsFactors=FALSE)
T6df$rse_T6 <- sigma
T6df$se.fit_T6 <- PGAM$se.fit

sites_to_fill <- unique(c(sites_of_area_bust, DDo$site_no[DDo$ed_rch_zone == 1]))
i <- 0
for(site in sites_to_fill) {
  for(decade in as.character(DDo$decade[DDo$site_no == site])) {
    i <- i + 1
    message(site," ", decade, " ", i)
    tmp <- data.frame(comid=DDo$comid[DDo$site_no == site & DDo$decade == decade],
                      site_no=site,
                      huc12=DDo$huc12[DDo$site_no == site & DDo$decade == decade],
                      decade=decade,
                      dec_long_va=DDo$dec_long_va[DDo$site_no == site & DDo$decade == decade],
                      dec_lat_va=DDo$dec_lat_va[DDo$site_no == site & DDo$decade == decade],
                      in_model_T6=FALSE,
                      T6=DDo$T6[DDo$site_no == site & DDo$decade == decade],
                      est_lwr_T6=NA, est_T6=NA, est_upr_T6=NA,
                      rse_T6=NA, se.fit_T6=NA, stringsAsFactors=FALSE)
    T6df <- rbind(T6df, tmp)
  }
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(T6, newdata=tmp, se.fit=TRUE)
}
T6df <- T6df[order(T6df$site_no, T6df$decade),]

for(site in sites_to_fill) {
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(T6, newdata=tmp, se.fit=TRUE)
  pgk <- gamIntervals(jnk, gam=T6, interval="prediction")
  df <- data.frame(comid=tmp$comid, site_no=tmp$site_no, huc12=tmp$huc12,
                   decade=tmp$decade,
                   dec_long_va=tmp$dec_long_va, dec_lat_va=tmp$dec_lat_va,
                   in_model_T6=0,
                   T6=tmp$T6,
                   est_lwr_T6=pgk$lwr,
                   est_T6    =pgk$fit,
                   est_upr_T6=pgk$upr, stringsAsFactors=FALSE)
  df$rse_T6 <- sigma
  df$se.fit_T6 <- pgk$se.fit
  T6df[T6df$site_no == site,] <- df
}

#write_feather(T6df, "all_gage_est_T6.feather")



save(bnd, DDo, DD, D, PPLO, L1, T2, T3, T4, T5, T6, file="Models.RData")

#source("fdcest_quantiles.R")

plot(PPLOdf$pplo, PPLOdf$est_pplo,
     xlab="Observed value", ylab="Fitted value")
mtext(paste0("PPLO: NSE=",round(NSE(PPLOdf$est_pplo, PPLOdf$pplo), digits=2)))
abline(0,1)

plot(log10(L1df$L1), log10(L1df$est_L1),
     xlab="Observed value", ylab="Fitted value")
mtext(paste0("L1: NSE=",round(NSE(L1df$bias_corr*L1df$est_L1, L1df$L1), digits=2)))
abline(0,1)

plot(T2df$T2, T2df$est_T2,
     xlab="Observed value", ylab="Fitted value")
mtext(paste0("T2: NSE=",round(NSE(T2df$est_T2, T2df$T2), digits=2)))
abline(0,1)

plot(T3df$T3, T3df$est_T3,
     xlab="Observed value", ylab="Fitted value")
mtext(paste0("T3: NSE=",round(NSE(T3df$est_T3, T3df$T3), digits=2)))
abline(0,1)

plot(T4df$T4, T4df$est_T4,
     xlab="Observed value", ylab="Fitted value")
mtext(paste0("T4: NSE=",round(NSE(T4df$est_T4, T4df$T4), digits=2)))
abline(0,1)

plot(T5df$T5, T5df$est_T5,
     xlab="Observed value", ylab="Fitted value")
mtext(paste0("T5: NSE=",round(NSE(T5df$est_T5, T5df$T5), digits=2)))
abline(0,1)

plot(T6df$T6, T6df$est_T6,
     xlab="Observed value", ylab="Fitted value")
mtext(paste0("T6: NSE=",round(NSE(T6df$est_T6, T6df$T6), digits=2)))
abline(0,1)

