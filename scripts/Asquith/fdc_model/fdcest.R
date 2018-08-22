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


"gamIntervals" <-
function(gam_predicts_with_se.fit, gam=NULL, sigma=NULL,
         interval=c("none", "confidence", "prediction"), level=0.95, ...) {
   # Demo: library(mgcv)
   #       X <- 2*pi*(1:360)/360 # simulate some X
   #       Y <- 1.6*sin(X) + 40*cos(X) + rnorm(length(X), sd=12)
   #     GAM <- gam(Y~s(X)); PGAM <- predict(GAM, se.fit=TRUE)
   #     PGAM <- gamIntervals(PGAM, gam=GAM, interval="confidence")
   #     print(head(PGAM))
   #     print(head(PGAM$leverage)); print(head(GAM$hat)) # see they are the value
   # plot(GAM$hat, (PGAM$se.fit/PGAM$residual.scale)^2)
   # Compare what the GAM says its leverage values are to back computed.
   # The plot() only works because predict() called back on the actual model.
   if(class(gam)[1] != "gam") {
      warning("need the actual GAM model too via the 'gam' argument")
      return()
   }
   z <- as.data.frame(gam_predicts_with_se.fit)
   if(! any(names(z) == "se.fit")) {
      warning("need gam predictions with se.fit=TRUE passed for 'gam_predicts_with_se.fit'")
      return()
   }
   interval <- match.arg(interval)
   sum.gam <- summary(gam); n <- sum.gam$n # summary.gam() and the sample size
   if(is.null(sigma)) sigma <- sqrt(sum.gam$scale)
   z$residual.scale <- sigma # residual standard error
   df <- n-sum(gam$edf)           # total degrees of freedom
   QT <- abs(qt((1-level)/2, df)) # will do the +/- separately
   z$leverage <- (z$se.fit/sigma)^2
   if(interval == "none") {
      z$lwr <- z$upr <- NA
   } else {
        one <- ifelse(interval == "confidence", 0, 1)
        tmp <- sqrt(one+z$leverage)
      z$lwr <- z$fit - sigma*QT*tmp
      z$upr <- z$fit + sigma*QT*tmp
   }
   attr(z, "interval")                  <- interval
   attr(z, "level")                     <- level
   attr(z, "t-dist_degrees_of_freedom") <- df
   return(z)
}


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
write.table(sitefile, file="sitefile.txt", quote=FALSE, row.names=FALSE, sep="\t")

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



DD$ppt_mean    <- log10(DD$ppt_mean)
DD$temp_mean   <- log10(DD$temp_mean)
DD$basin_area  <- log10(DD$basin_area)
DD$basin_slope <- log10(DD$basin_slope/100)

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


plot(DD$CDA, DD$basin_area, lwd=0.5,
     xlab="log10(NWIS CDA)", ylab="log10(NHDplus basin area")
abline(0,1)
abline(1/3,1,lty=2); abline(-1/3,1,lty=2)
abline(1/2,1,lty=2); abline(-1/2,1,lty=3)
mtext("Diagnostic check on watershed areas")
jnk <- abs(DD$basin_area - DD$CDA)
summary(jnk[jnk > 1/2])
sites_of_area_bust <- unique(DD$site_no[jnk > 1/2])
DD_sites_of_area_bust <- DD[DD$site_no == sites_of_area_bust[1], ]
for(site in sites_of_area_bust[2:length(sites_of_area_bust)]) {
  DD_sites_of_area_bust <- rbind(DD_sites_of_area_bust,DD[DD$site_no == site, ] )
}
for(site in sites_of_area_bust) {
  points(DD$CDA[DD$site_no == site], DD$basin_area[DD$site_no == site], pch=16, col=2)
  DD <- DD[DD$site_no != site,]
}
text(0,5, paste(sites_of_area_bust, collapse=", "), cex=0.6, pos=4)
#points(DD$CDA[DD$site_no == "08167000"],
#       DD$basin_area[DD$site_no == "08167000"], pch=16, col=4)




D <- DD;

D <- D[D$ed_rch_zone != "1",]

duan_smearing_estimator <- function(model) { sum(10^residuals(model))/length(residuals(model)) }

save(bnd, D, DD, DDo, knots, bnd,
     DD_sites_of_area_bust, duan_smearing_estimator, file="DEMO.RData")

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
#[1] 2007
length(Z$n[Z$nzero != 0])
#[1] 739
length(Pto[Pto < log10(3653)])
#[1] 387
length(Gto[Gto < log10(3653)])
#[1] 414
txt <- paste0("Number of streamgage:decade records with no zero-flow conditions: 2,011\n",
              "Number of streamgage:decade records with zero-flow conditions: 738\n",
              "Number of survival regression predicted with zero-flow conditions: 388\n",
              "Number of GAM regression predicted with zero-flow conditions: 409\n")
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

GM1 <- gam(flowtime~basin_area+
                    s(ppt_mean, k=5)+s(temp_mean, k=4)+s(dni_ann, k=7)+
                    developed+s(grassland)+
                    bedperm+decade-1,
           family=tobit1(left.threshold=  Z$left.threshold,
                         right.threshold=Z$right.threshold), data=Z); summary(GM1)
GM2 <- gam(flowtime~basin_area+
                    s(ppt_mean, k=5)+s(temp_mean, k=4)+s(dni_ann, k=7)+
                    developed+s(grassland)+
                    bedperm+decade-1+
        s(x,y, bs="so", xt=list(bnd=bnd)), knots=knots,
           family=tobit1(left.threshold=  Z$left.threshold,
                         right.threshold=Z$right.threshold), data=Z); summary(GM2)
GM3 <- gam(flowtime~basin_area+
                    s(ppt_mean, k=5)+s(temp_mean, k=4)+s(dni_ann, k=7)+
                    developed+s(grassland)+
                    bedperm+decade-1+
        s(x,y), knots=knots,
           family=tobit1(left.threshold=  Z$left.threshold,
                         right.threshold=Z$right.threshold), data=Z); summary(GM3)

P <- Po <- predict(SM0); P[P > log10(3653)] <- log10(3653)
C1 <- C1o <- predict(GM1); C1[C1 > log10(3653)] <- log10(3653)
C2 <- C2o <- predict(GM2); C2[C2 > log10(3653)] <- log10(3653)
C3 <- C3o <- predict(GM3); C3[C3 > log10(3653)] <- log10(3653)
pa <- abs(P - Z$flowtime)
ca <- abs(C1- Z$flowtime)
cb <- abs(C2- Z$flowtime)
cc <- abs(C3- Z$flowtime)
summary(pa)
summary(ca)
summary(cb)
summary(cc)


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
points(qnorm(pp(cc)), sort(cc), col=3, lwd=0.5)

Ppplo  <- (3653-10^P) /3653
C1pplo <- (3653-10^C1)/3653
C2pplo <- (3653-10^C2)/3653
C3pplo <- (3653-10^C3)/3653
plot(  Z$pplo, Ppplo, xlim=c(0,1), ylim=c(0,1))
points(Z$pplo, C2pplo, col=2)
points(Z$pplo, C2pplo, col=4)
points(Z$pplo, C3pplo, col=3)
abline(0,1)

PPLO <- GM3
save(DDo, DD, D, SM0, GM1, GM2, PPLO, file="PPLOS.RData")

pplo.sigma <- sqrt(mean((predict(PPLO) - Z$flowtime)^2))
PGAM <- gamIntervals(predict(GM3, se.fit=TRUE), gam=GM3, interval="prediction", sigma=pplo.sigma)
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
  pgk <- gamIntervals(jnk, gam=GM3, interval="prediction", sigma=pplo.sigma)
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
PPLOdf$in_model_pplo[PPLOdf$in_model_pplo == 1] <- "yes"
PPLOdf$in_model_pplo[PPLOdf$in_model_pplo == 0] <- "no"
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

write_feather(PPLOdf, "all_gage_est_pplo.feather")



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


###### END PPLO


Z <- D
x <- Z$x; y <- Z$y
z <- log10(Z$L1)
L1   <- gam(z~basin_area +
              s(ppt_mean, k=5) + s(temp_mean, k=4) + s(dni_ann, k=7)+
              developed+
              mixed_forest+shrubland+
              decade-1+
              s(x, y), knots=knots, data=Z, # , bs="so", xt=list(bnd=bnd)
              family="gaussian")
L1$duan_smearing <- duan_smearing_estimator(L1)
pdf("L1.pdf", useDingbats=FALSE)
  plot(z, fitted.values(L1))
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
# Terms invert in upper/lower meaning and hence the flipping during data.frame construction.
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

write_feather(L1df, "all_gage_est_L1.feather")



z <- D$L2/D$L1 # --------------------------- Coefficient of L-variation
T2   <- gam(z~basin_area +
              s(ppt_mean, k=5) + s(temp_mean, k=4) + s(dni_ann, k=7)+
              developed+
              mixed_forest+shrubland+
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

write_feather(T2df, "all_gage_est_T2.feather")




z <- D$T3      # --------------------------- L-skew
T3   <- gam(z~basin_area +
              s(ppt_mean, k=5) + s(temp_mean, k=4) + s(dni_ann, k=7)+
              developed+
              mixed_forest+shrubland+
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

write_feather(T3df, "all_gage_est_T3.feather")





z <- D$T4      # --------------------------- L-kurtosis
T4   <- gam(z~basin_area +
              s(ppt_mean, k=5) + s(temp_mean, k=4) + s(dni_ann, k=7)+
              developed+
              mixed_forest+shrubland+
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

write_feather(T4df, "all_gage_est_T4.feather")




z <- D$T5      # --------------------------- Fifth L-moment ratio
T5   <- gam(z~basin_area +
              s(ppt_mean, k=5) + s(temp_mean, k=4) + s(dni_ann, k=7)+
              developed+
              mixed_forest+shrubland+
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

write_feather(T5df, "all_gage_est_T5.feather")





z <- D$T6      # --------------------------- Sixth L-moment ratio
T6   <- gam(z~basin_area +
              s(ppt_mean, k=5) + s(temp_mean, k=4) + s(dni_ann, k=7)+
              developed+
              mixed_forest+shrubland+
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

write_feather(T6df, "all_gage_est_T6.feather")






z <- log10(D$f50+1)      # ---------------------------
F50   <- gam(z~basin_area +
               s(ppt_mean, k=5) + s(temp_mean, k=4) + s(dni_ann, k=7)+
               developed+
               mixed_forest+shrubland+
               decade-1+
               s(x, y), knots=knots, data=D, #, bs="so", xt=list(bnd=bnd)
             family="gaussian")
pdf("F50.pdf", useDingbats=FALSE)
plot(z, fitted.values(F50))
abline(0,1)
plot(F50, scheme=2, residuals=TRUE)
points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
points(knots$x, knots$y, pch=16, cex=1.1, col=4)
text(100, 500, "F50")
dev.off()



PGAM <- gamIntervals(predict(F50, se.fit=TRUE), gam=F50, interval="prediction")
#plot(F50$hat, (PGAM$se.fit/PGAM$residual.scale)^2)
#abline(0,1)
sigma <- PGAM$residual.scale[1]
# Terms invert in upper/lower meaning and hence the flipping during data.frame construction.
F50df <- data.frame(comid=Z$comid, site_no=Z$site_no, huc12=Z$huc12,
                     decade=Z$decade,
                     dec_long_va=Z$dec_long_va, dec_lat_va=Z$dec_lat_va,
                     in_model_f50=TRUE, f50=Z$f50,
                     est_lwr_f50=10^PGAM$lwr-1,
                     est_f50    =10^PGAM$fit-1,
                     est_upr_f50=10^PGAM$upr-1, stringsAsFactors=FALSE)
F50df$rse_f50 <- sigma
F50df$se.fit_f50 <- PGAM$se.fit

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
                      in_model_f50=FALSE,
                      f50=DDo$f50[DDo$site_no == site & DDo$decade == decade],
                      est_lwr_f50=NA, est_f50=NA, est_upr_f50=NA,
                      rse_f50=NA, se.fit_f50=NA, stringsAsFactors=FALSE)
    F50df <- rbind(F50df, tmp)
  }
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(F50, newdata=tmp, se.fit=TRUE)
}
F50df <- F50df[order(F50df$site_no, F50df$decade),]

for(site in sites_to_fill) {
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(F50, newdata=tmp, se.fit=TRUE)
  pgk <- gamIntervals(jnk, gam=F50, interval="prediction")
  df <- data.frame(comid=tmp$comid, site_no=tmp$site_no, huc12=tmp$huc12,
                   decade=tmp$decade,
                   dec_long_va=tmp$dec_long_va, dec_lat_va=tmp$dec_lat_va,
                   in_model_f50=0,
                   f50=tmp$f50,
                   est_lwr_f50=10^pgk$lwr-1,
                   est_f50    =10^pgk$fit-1,
                   est_upr_f50=10^pgk$upr-1, stringsAsFactors=FALSE)
  df$rse_f50 <- sigma
  df$se.fit_f50 <- pgk$se.fit
  F50df[F50df$site_no == site,] <- df
}
F50df$est_f50[F50df$est_f50 < 0] <- 0
F50df$est_lwr_f50[F50df$est_lwr_f50 < 0] <- 0
F50df$est_upr_f50[F50df$est_upr_f50 < 0] <- 0

write_feather(F50df, "all_gage_est_f50.feather")



#sink("right_tail_flowing_fdc.txt")
z <- log10(D$f90+1)      # ---------------------------
F90   <- gam(z~basin_area  + s(flood_storage, k=7) +
              s(ppt_mean, k=5) + s(dni_ann, k=7)+
              developed+
              mixed_forest+shrubland+
              decade-1+
              s(x, y), knots=knots, data=D, #, bs="so", xt=list(bnd=bnd)
              family="gaussian")
print(summary(F90))
pdf("F90.pdf", useDingbats=FALSE)
  plot(z, fitted.values(F90))
  abline(0,1)
  plot(F90, scheme=2, residuals=TRUE)
  points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
  points(knots$x, knots$y, pch=16, cex=1.1, col=4)
  text(100, 500, "F90")
dev.off()



PGAM <- gamIntervals(predict(F90, se.fit=TRUE), gam=F90, interval="prediction")
#plot(F90$hat, (PGAM$se.fit/PGAM$residual.scale)^2)
#abline(0,1)
sigma <- PGAM$residual.scale[1]
# Terms invert in upper/lower meaning and hence the flipping during data.frame construction.
F90df <- data.frame(comid=Z$comid, site_no=Z$site_no, huc12=Z$huc12,
                     decade=Z$decade,
                     dec_long_va=Z$dec_long_va, dec_lat_va=Z$dec_lat_va,
                     in_model_f90=TRUE, f90=Z$f90,
                     est_lwr_f90=10^PGAM$lwr-1,
                     est_f90    =10^PGAM$fit-1,
                     est_upr_f90=10^PGAM$upr-1, stringsAsFactors=FALSE)
F90df$rse_f90 <- sigma
F90df$se.fit_f90 <- PGAM$se.fit

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
                      in_model_f90=FALSE,
                      f90=DDo$f90[DDo$site_no == site & DDo$decade == decade],
                      est_lwr_f90=NA, est_f90=NA, est_upr_f90=NA,
                      rse_f90=NA, se.fit_f90=NA, stringsAsFactors=FALSE)
    F90df <- rbind(F90df, tmp)
  }
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(F90, newdata=tmp, se.fit=TRUE)
}
F90df <- F90df[order(F90df$site_no, F90df$decade),]

for(site in sites_to_fill) {
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(F90, newdata=tmp, se.fit=TRUE)
  pgk <- gamIntervals(jnk, gam=F90, interval="prediction")
  df <- data.frame(comid=tmp$comid, site_no=tmp$site_no, huc12=tmp$huc12,
                   decade=tmp$decade,
                   dec_long_va=tmp$dec_long_va, dec_lat_va=tmp$dec_lat_va,
                   in_model_f90=0,
                   f90=tmp$f90,
                   est_lwr_f90=10^pgk$lwr-1,
                   est_f90    =10^pgk$fit-1,
                   est_upr_f90=10^pgk$upr-1, stringsAsFactors=FALSE)
  df$rse_f90 <- sigma
  df$se.fit_f90 <- pgk$se.fit
  F90df[F90df$site_no == site,] <- df
}
F90df$est_f90[F90df$est_f90 < 0] <- 0
F90df$est_lwr_f90[F90df$est_lwr_f90 < 0] <- 0
F90df$est_upr_f90[F90df$est_upr_f90 < 0] <- 0


write_feather(F90df, "all_gage_est_f90.feather")



z <- log10(D$f95+1)      # ---------------------------
F95   <- gam(z~basin_area + s(flood_storage, k=7) +
               s(ppt_mean, k=5) + s(dni_ann, k=7)+
               developed+
               mixed_forest+shrubland+
               decade-1+
               s(x, y), knots=knots, data=D, #, bs="so", xt=list(bnd=bnd)
             family="gaussian")
print(summary(F95))
pdf("F95.pdf", useDingbats=FALSE)
plot(z, fitted.values(F95))
abline(0,1)
plot(F95, scheme=2, residuals=TRUE)
points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
points(knots$x, knots$y, pch=16, cex=1.1, col=4)
text(100, 500, "F95")
dev.off()



PGAM <- gamIntervals(predict(F95, se.fit=TRUE), gam=F95, interval="prediction")
#plot(F95$hat, (PGAM$se.fit/PGAM$residual.scale)^2)
#abline(0,1)
sigma <- PGAM$residual.scale[1]
# Terms invert in upper/lower meaning and hence the flipping during data.frame construction.
F95df <- data.frame(comid=Z$comid, site_no=Z$site_no, huc12=Z$huc12,
                     decade=Z$decade,
                     dec_long_va=Z$dec_long_va, dec_lat_va=Z$dec_lat_va,
                     in_model_f95=TRUE, f95=Z$f95,
                     est_lwr_f95=10^PGAM$lwr-1,
                     est_f95    =10^PGAM$fit-1,
                     est_upr_f95=10^PGAM$upr-1, stringsAsFactors=FALSE)
F95df$rse_f95 <- sigma
F95df$se.fit_f95 <- PGAM$se.fit

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
                      in_model_f95=FALSE,
                      f95=DDo$f95[DDo$site_no == site & DDo$decade == decade],
                      est_lwr_f95=NA, est_f95=NA, est_upr_f95=NA,
                      rse_f95=NA, se.fit_f95=NA, stringsAsFactors=FALSE)
    F95df <- rbind(F95df, tmp)
  }
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(F95, newdata=tmp, se.fit=TRUE)
}
F95df <- F95df[order(F95df$site_no, F95df$decade),]

for(site in sites_to_fill) {
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(F95, newdata=tmp, se.fit=TRUE)
  pgk <- gamIntervals(jnk, gam=F95, interval="prediction")
  df <- data.frame(comid=tmp$comid, site_no=tmp$site_no, huc12=tmp$huc12,
                   decade=tmp$decade,
                   dec_long_va=tmp$dec_long_va, dec_lat_va=tmp$dec_lat_va,
                   in_model_f95=0,
                   f95=tmp$f95,
                   est_lwr_f95=10^pgk$lwr-1,
                   est_f95    =10^pgk$fit-1,
                   est_upr_f95=10^pgk$upr-1, stringsAsFactors=FALSE)
  df$rse_f95 <- sigma
  df$se.fit_f95 <- pgk$se.fit
  F95df[F95df$site_no == site,] <- df
}
F95df$est_f95[F95df$est_f95 < 0] <- 0
F95df$est_lwr_f95[F95df$est_lwr_f95 < 0] <- 0
F95df$est_upr_f95[F95df$est_upr_f95 < 0] <- 0

write_feather(F95df, "all_gage_est_f95.feather")




z <- log10(D$f98+1)      # ---------------------------
F98   <- gam(z~basin_area + s(flood_storage, k=7) +
               s(ppt_mean, k=5) + s(dni_ann, k=7)+
               developed+
               mixed_forest+shrubland+
               decade-1+
               s(x, y), knots=knots, data=D, #, bs="so", xt=list(bnd=bnd)
             family="gaussian")
print(summary(F98))
pdf("F98.pdf", useDingbats=FALSE)
plot(z, fitted.values(F98))
abline(0,1)
plot(F98, scheme=2, residuals=TRUE)
points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
points(knots$x, knots$y, pch=16, cex=1.1, col=4)
text(100, 500, "F98")
dev.off()

PGAM <- gamIntervals(predict(F98, se.fit=TRUE), gam=F98, interval="prediction")
#plot(F98$hat, (PGAM$se.fit/PGAM$residual.scale)^2)
#abline(0,1)
sigma <- PGAM$residual.scale[1]
# Terms invert in upper/lower meaning and hence the flipping during data.frame construction.
F98df <- data.frame(comid=Z$comid, site_no=Z$site_no, huc12=Z$huc12,
                     decade=Z$decade,
                     dec_long_va=Z$dec_long_va, dec_lat_va=Z$dec_lat_va,
                     in_model_f98=TRUE, f98=Z$f98,
                     est_lwr_f98=10^PGAM$lwr-1,
                     est_f98    =10^PGAM$fit-1,
                     est_upr_f98=10^PGAM$upr-1, stringsAsFactors=FALSE)
F98df$rse_f98 <- sigma
F98df$se.fit_f98 <- PGAM$se.fit

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
                      in_model_f98=FALSE,
                      f98=DDo$f98[DDo$site_no == site & DDo$decade == decade],
                      est_lwr_f98=NA, est_f98=NA, est_upr_f98=NA,
                      rse_f98=NA, se.fit_f98=NA, stringsAsFactors=FALSE)
    F98df <- rbind(F98df, tmp)
  }
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(F98, newdata=tmp, se.fit=TRUE)
}
F98df <- F98df[order(F98df$site_no, F98df$decade),]

for(site in sites_to_fill) {
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(F98, newdata=tmp, se.fit=TRUE)
  pgk <- gamIntervals(jnk, gam=F98, interval="prediction")
  df <- data.frame(comid=tmp$comid, site_no=tmp$site_no, huc12=tmp$huc12,
                   decade=tmp$decade,
                   dec_long_va=tmp$dec_long_va, dec_lat_va=tmp$dec_lat_va,
                   in_model_f98=0,
                   f98=tmp$f98,
                   est_lwr_f98=10^pgk$lwr-1,
                   est_f98    =10^pgk$fit-1,
                   est_upr_f98=10^pgk$upr-1, stringsAsFactors=FALSE)
  df$rse_f98 <- sigma
  df$se.fit_f98 <- pgk$se.fit
  F98df[F98df$site_no == site,] <- df
}
F98df$est_f98[F98df$est_f98 < 0] <- 0
F98df$est_lwr_f98[F98df$est_lwr_f98 < 0] <- 0
F98df$est_upr_f98[F98df$est_upr_f98 < 0] <- 0

write_feather(F98df, "all_gage_est_f98.feather")




z <- log10(D$f99+1)      # ---------------------------
F99   <- gam(z~basin_area + s(flood_storage, k=7) +
               s(ppt_mean, k=5) + s(dni_ann, k=7)+
               developed+
               mixed_forest+shrubland+
               decade-1+
               s(x, y), knots=knots, data=D, #, bs="so", xt=list(bnd=bnd)
             family="gaussian")
print(summary(F99))
pdf("F99.pdf", useDingbats=FALSE)
plot(z, fitted.values(F99))
abline(0,1)
plot(F99, scheme=2, residuals=TRUE)
points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
points(knots$x, knots$y, pch=16, cex=1.1, col=4)
text(100, 500, "F99")
dev.off()


PGAM <- gamIntervals(predict(F99, se.fit=TRUE), gam=F99, interval="prediction")
#plot(F99$hat, (PGAM$se.fit/PGAM$residual.scale)^2)
#abline(0,1)
sigma <- PGAM$residual.scale[1]
# Terms invert in upper/lower meaning and hence the flipping during data.frame construction.
F99df <- data.frame(comid=Z$comid, site_no=Z$site_no, huc12=Z$huc12,
                     decade=Z$decade,
                     dec_long_va=Z$dec_long_va, dec_lat_va=Z$dec_lat_va,
                     in_model_f99=TRUE, f99=Z$f99,
                     est_lwr_f99=10^PGAM$lwr-1,
                     est_f99    =10^PGAM$fit-1,
                     est_upr_f99=10^PGAM$upr-1, stringsAsFactors=FALSE)
F99df$rse_f99 <- sigma
F99df$se.fit_f99 <- PGAM$se.fit

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
                      in_model_f99=FALSE,
                      f99=DDo$f99[DDo$site_no == site & DDo$decade == decade],
                      est_lwr_f99=NA, est_f99=NA, est_upr_f99=NA,
                      rse_f99=NA, se.fit_f99=NA, stringsAsFactors=FALSE)
    F99df <- rbind(F99df, tmp)
  }
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(F99, newdata=tmp, se.fit=TRUE)
}
F99df <- F99df[order(F99df$site_no, F99df$decade),]

for(site in sites_to_fill) {
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(F99, newdata=tmp, se.fit=TRUE)
  pgk <- gamIntervals(jnk, gam=F99, interval="prediction")
  df <- data.frame(comid=tmp$comid, site_no=tmp$site_no, huc12=tmp$huc12,
                   decade=tmp$decade,
                   dec_long_va=tmp$dec_long_va, dec_lat_va=tmp$dec_lat_va,
                   in_model_f99=0,
                   f99=tmp$f99,
                   est_lwr_f99=10^pgk$lwr-1,
                   est_f99    =10^pgk$fit-1,
                   est_upr_f99=10^pgk$upr-1, stringsAsFactors=FALSE)
  df$rse_f99 <- sigma
  df$se.fit_f99 <- pgk$se.fit
  F99df[F99df$site_no == site,] <- df
}
F99df$est_f99[F99df$est_f99 < 0] <- 0
F99df$est_lwr_f99[F99df$est_lwr_f99 < 0] <- 0
F99df$est_upr_f99[F99df$est_upr_f99 < 0] <- 0

write_feather(F99df, "all_gage_est_f99.feather")



z <- log10(D$f99.9+1)      # ---------------------------
F99p9   <- gam(z~basin_area + s(flood_storage, k=7) +
               s(ppt_mean, k=5) + s(dni_ann, k=7)+
               developed+
               mixed_forest+shrubland+
               decade-1+
               s(x, y), knots=knots, data=D, #, bs="so", xt=list(bnd=bnd)
             family="gaussian")
print(summary(F99p9))
pdf("F99p9.pdf", useDingbats=FALSE)
plot(z, fitted.values(F99p9))
abline(0,1)
plot(F99p9, scheme=2, residuals=TRUE)
points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
points(knots$x, knots$y, pch=16, cex=1.1, col=4)
text(100, 500, "F99p9")
dev.off()



PGAM <- gamIntervals(predict(F99p9, se.fit=TRUE), gam=F99p9, interval="prediction")
#plot(F99p9$hat, (PGAM$se.fit/PGAM$residual.scale)^2)
#abline(0,1)
sigma <- PGAM$residual.scale[1]
# Terms invert in upper/lower meaning and hence the flipping during data.frame construction.
F99p9df <- data.frame(comid=Z$comid, site_no=Z$site_no, huc12=Z$huc12,
                     decade=Z$decade,
                     dec_long_va=Z$dec_long_va, dec_lat_va=Z$dec_lat_va,
                     in_model_f99.9=TRUE, f99.9=Z$f99.9,
                     est_lwr_f99.9=10^PGAM$lwr-1,
                     est_f99.9    =10^PGAM$fit-1,
                     est_upr_f99.9=10^PGAM$upr-1, stringsAsFactors=FALSE)
F99p9df$rse_f99.9 <- sigma
F99p9df$se.fit_f99.9 <- PGAM$se.fit

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
                      in_model_f99.9=FALSE,
                      f99.9=DDo$f99.9[DDo$site_no == site & DDo$decade == decade],
                      est_lwr_f99.9=NA, est_f99.9=NA, est_upr_f99.9=NA,
                      rse_f99.9=NA, se.fit_f99.9=NA, stringsAsFactors=FALSE)
    F99p9df <- rbind(F99p9df, tmp)
  }
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(F99p9, newdata=tmp, se.fit=TRUE)
}
F99p9df <- F99p9df[order(F99p9df$site_no, F99p9df$decade),]

for(site in sites_to_fill) {
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(F99p9, newdata=tmp, se.fit=TRUE)
  pgk <- gamIntervals(jnk, gam=F99p9, interval="prediction")
  df <- data.frame(comid=tmp$comid, site_no=tmp$site_no, huc12=tmp$huc12,
                   decade=tmp$decade,
                   dec_long_va=tmp$dec_long_va, dec_lat_va=tmp$dec_lat_va,
                   in_model_f99.9=0,
                   f99.9=tmp$f99.9,
                   est_lwr_f99.9=10^pgk$lwr-1,
                   est_f99.9    =10^pgk$fit-1,
                   est_upr_f99.9=10^pgk$upr-1, stringsAsFactors=FALSE)
  df$rse_f99.9 <- sigma
  df$se.fit_f99.9 <- pgk$se.fit
  F99p9df[F99p9df$site_no == site,] <- df
}
F99p9df$est_f99.9[F99p9df$est_f99.9 < 0] <- 0
F99p9df$est_lwr_f99.9[F99p9df$est_lwr_f99.9 < 0] <- 0
F99p9df$est_upr_f99.9[F99p9df$est_upr_f99.9 < 0] <- 0

write_feather(F99p9df, "all_gage_est_f99p9.feather")

#sink()

save(DDo, DD, D, PPLO, L1, T2, T3, T4, T5, T6,
     F50, F90, F95, F98, F99, F99p9, file="Models.RData")


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

plot(log10(F50df$f50), log10(F50df$est_f50),
     xlab="Observed value", ylab="Fitted value")
mtext(paste0("f50: NSE=",round(NSE(F50df$est_f50, F50df$f50), digits=2)))
abline(0,1)

plot(log10(F90df$f90), log10(F90df$est_f90),
     xlab="Observed value", ylab="Fitted value")
mtext(paste0("f90: NSE=",round(NSE(F90df$est_f90, F90df$f90), digits=2)))
abline(0,1)

plot(log10(F95df$f95), log10(F95df$est_f95),
     xlab="Observed value", ylab="Fitted value")
mtext(paste0("f95: NSE=",round(NSE(F95df$est_f95, F95df$f95), digits=2)))
abline(0,1)

plot(log10(F98df$f98), log10(F98df$est_f98),
     xlab="Observed value", ylab="Fitted value")
mtext(paste0("f98: NSE=",round(NSE(F98df$est_f98, F98df$f98), digits=2)))
abline(0,1)

plot(log10(F99df$f99), log10(F99df$est_f99),
     xlab="Observed value", ylab="Fitted value")
mtext(paste0("f99: NSE=",round(NSE(F99df$est_f99, F99df$f99), digits=2)))
abline(0,1)

plot(log10(F99p9df$f99.9), log10(F99p9df$est_f99.9),
     xlab="Observed value", ylab="Fitted value")
mtext(paste0("f99.9: NSE=",round(NSE(F99p9df$est_f99.9, F99p9df$f99.9), digits=2)))
abline(0,1)

#"F50df"   "F90df"   "F95df"   "F98df"   "F99df"   "F99p9df" "L1df"    "PPLOdf"
#"T2df"    "T3df"    "T4df"    "T5df"    "T6df"


cvPPLOgo <- function(parent=PPLOdf, sigma=pplo.sigma, sites_to_fill=sites_to_fill) {
  sites <- decades <- values <- ests <- estlwrs <- estuprs <- NULL
  estft <- estlwrft <- estuprft <- NULL
  rse <- rsecv <- biascv <- 0; i <- 0
  for(site in unique(D$site_no)) { i <- i + 1
    Z <- D[D$site_no != site,]; x <- Z$x; y <- Z$y
    Z$left.threshold <-  log10(rep(0, length(Z$nzero)))
    Z$right.threshold <- log10(Z$n)
    Z$flowtime <- log10(Z$n - Z$nzero)
    model <- gam(flowtime~basin_area+
                    s(ppt_mean, k=5)+s(temp_mean, k=4)+s(dni_ann, k=7)+
                    developed+s(grassland)+
                    bedperm+decade-1+
        s(x,y), knots=knots,
           family=tobit1(left.threshold=  Z$left.threshold,
                         right.threshold=Z$right.threshold), data=Z)
    rse[i] <- sigma
    new.sigma <- sqrt(mean((predict(model)-Z$flowtime)^2))
    tmp <- D[D$site_no == site,]; val <- tmp$pplo # LOOK HERE
    pgk <- predict(model, newdata=tmp, se.fit=TRUE)
    pgk <- gamIntervals(pgk, gam=model, interval="prediction", sigma=new.sigma)
    pp  <- (3653-10^pgk$fit)/3653; pp[ pp  < 0] <- 0
    ppl <- (3653-10^pgk$lwr)/3653; ppl[ppl < 0] <- 0
    ppu <- (3653-10^pgk$upr)/3653; ppu[ppu < 0] <- 0
    ft  <- pgk$fit; lft <- pgk$lwr; uft <- pgk$upr
    res <- val - pp; biascv[i] <- mean(res)
    res <- sum(res^2); rsecv[i] <- new.sigma
    message(site, ", Bias=",biascv[i],", RSE=",rse[i],", RSEcv=",rsecv[i])
    if(is.null(values)) {
      sites <- tmp$site_no
      decades <- tmp$decade
      values <- val
      ests <- pp
      estlwrs <- ppu # swap is correct
      estuprs <- ppl # swap is correct
      estft <- ft; estlwrft <- lft; estuprft <- uft
    } else {
      s <- (length(values)+1); ss <- s:(s+length(val)-1)
      sites[ss] <- tmp$site_no
      decades[ss] <- tmp$decade
      values[ss] <- val
      ests[ss] <- pp
      estlwrs[ss] <- ppu # swap is correct
      estuprs[ss] <- ppl # swap is correct
      estft[ss] <- ft; estlwrft[ss] <- lft; estuprft[ss] <- uft
    }
    #print(parent$est_lwr_pplo[parent$site_no == site])
    #print(val)
    #print(parent$est_upr_pplo[parent$site_no == site])
  }
  print(summary(biascv));   print(summary(abs(biascv)))
  print(summary(rse))
  print(summary(rsecv))
  zz <- data.frame(site_no_bak=sites, decade_bak=decades, orgfit=values,
             loo_est_lwr_pplo=estlwrs, loo_est_pplo=ests, loo_est_upr_pplo=estuprs,
             loo_est_lwr_flowtime=estlwrft, loo_est_flowtime=estft, loo_est_upr_flowtime=estuprft,
             stringsAsFactors=FALSE)
  return(zz)
}

cvL1go <- function(parent=L1df) {
  sites <- decades <- values <- ests <- estlwrs <- estuprs <- NULL
  duans <- NULL
  rse <- rsecv <- biascv <- 0; i <- 0
  for(site in unique(D$site_no)) { i <- i + 1
    Z <- D[D$site_no != site,]; x <- Z$x; y <- Z$y
    z <- log10(Z$L1) # LOOK HERE
    model   <- gam(z~basin_area +
               s(ppt_mean, k=5) + s(temp_mean, k=4) + s(dni_ann, k=7)+
               developed+
               mixed_forest+shrubland+
               decade-1+
               s(x, y), knots=knots, data=Z, #, bs="so", xt=list(bnd=bnd)
               family="gaussian")
    duan <- duan_smearing_estimator(model); #print(duan)
    rse[i] <- sqrt(summary(model)$scale)
    tmp <- D[D$site_no == site,]; val <- log10(tmp$L1) # LOOK HERE
    pgk <- predict(model, newdata=tmp, se.fit=TRUE)
    pgk <- gamIntervals(pgk, gam=model, interval="prediction")
    res <- (val - pgk$fit); biascv[i] <- mean(res)
    res <- sum(res^2); rsecv[i] <- sqrt(mean(res))
    message(site, ", Bias=",biascv[i],", RSE=",rse[i],", RSEcv=",rsecv[i])
    if(is.null(values)) {
      sites <- tmp$site_no
      decades <- tmp$decade
      duans <- rep(duan, length(tmp$site_no))
      values <- val
      ests <- pgk$fit
      estlwrs <- pgk$lwr
      estuprs <- pgk$upr
    } else {
      s <- (length(values)+1); ss <- s:(s+length(val)-1)
      sites[ss] <- tmp$site_no
      decades[ss] <- tmp$decade
      duans[ss] <- rep(duan, length(val))
      values[ss] <- val
      ests[ss] <- pgk$fit
      estlwrs[ss] <- pgk$lwr
      estuprs[ss] <- pgk$upr
    }
    #print(parent$est_lwr_L1[parent$site_no == site])
    #print(10^val)
    #print(parent$est_upr_L1[parent$site_no == site])
  }
  print(summary(biascv));   print(summary(abs(biascv)))
  print(summary(rse))
  print(summary(rsecv))
  zz <- data.frame(site_no_bak=sites, decade_bak=decades, loo_bias_corr=duans, orgfit=values,
             loo_est_lwr_L1=10^estlwrs, loo_est_L1=10^ests, loo_est_upr_L1=10^estuprs,
             stringsAsFactors = FALSE)
  return(zz)
}

cvL1 <- cvL1go()
cvPPLO <- cvPPLOgo()

  for(site in sites_to_fill) {
     tmp <- DDo[DDo$site_no == site,]
     df <- data.frame(site_no_bak=tmp$site_no, decade_bak=tmp$decade,
                      loo_bias_corr=NA, orgfit=tmp$L1,
                      loo_est_lwr_L1=NA, loo_est_L1=NA, loo_est_upr_L1=NA,
                      stringsAsFactors=FALSE)
     cvL1 <- rbind(cvL1, df)
  }
  cvL1 <- cvL1[order(cvL1$site_no_bak, cvL1$decade_bak),]
  cvL1$key <-  paste(cvL1$site_no_bak, cvL1$decade_bak, sep=":")


L1df$key <- paste(L1df$site_no, L1df$decade, sep=":")
L1df_loo <- merge(L1df, cvL1, add=TRUE)
L1df_loo$key         <- NULL
L1df_loo$site_no_bak <- NULL
L1df_loo$decade_bak  <- NULL
L1df_loo$orgfit      <- NULL

  for(site in sites_to_fill) {
     tmp <- DDo[DDo$site_no == site,]
     df <- data.frame(site_no_bak=tmp$site_no, decade_bak=tmp$decade,
                      orgfit=tmp$pplo,
                      loo_est_lwr_pplo=NA, loo_est_pplo=NA, loo_est_upr_pplo=NA,
                      loo_est_lwr_flowtime=NA, loo_est_flowtime=NA, loo_est_upr_flowtime=NA,
                      stringsAsFactors=FALSE)
     cvPPLO <- rbind(cvPPLO, df)
  }
  cvPPLO <- cvPPLO[order(cvPPLO$site_no_bak, cvPPLO$decade_bak),]
  cvPPLO$key <-  paste(cvPPLO$site_no_bak, cvPPLO$decade_bak, sep=":")

PPLOdf$key <- paste(PPLOdf$site_no, PPLOdf$decade, sep=":")
PPLOdf_loo <- merge(PPLOdf, cvPPLO, add=TRUE)
PPLOdf_loo$key         <- NULL
PPLOdf_loo$site_no_bak <- NULL
PPLOdf_loo$decade_bak  <- NULL
PPLOdf_loo$orgfit      <- NULL


TrueMean <- 0*PPLOdf_loo$loo_est_pplo + (1-PPLOdf_loo$loo_est_pplo)*L1df_loo$loo_est_L1
TrueMean[L1df_loo$site_no == "08167000"]

