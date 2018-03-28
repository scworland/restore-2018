library(feather)
library(proxy)
library(akqdecay)
library(mgcv)
library(sp)
library(rgeos)
library(lmomco)
library(Lmoments)
library(survival)
# https://www.sciencebase.gov/catalog/item/5669a79ee4b08895842a1d47
FDC <- read_feather(file.choose()) # "all_gage_data.feather"
load(file.choose()) # "spRESTORE_MGCV_BND.RData"
load(file.choose()) # spDNI_1998to2009.RData

sites <- unique(FDC$site_no)
sitefile <- dataRetrieval::readNWISsite(sites)
sitefile <- sitefile[sitefile$agency_cd != "USCE",]

CDA <- sitefile$contrib_drain_area_va
CDA[is.na(CDA)] <- sitefile$drain_area_va[is.na(CDA)]
CDA <- pmin(sitefile$drain_area_va, sitefile$contrib_drain_area_va, na.rm=TRUE)
sitefile$contrib_drain_area_va <- CDA <- CDA*2.589988 # move to km2


LATLONG <- paste0("+init=epsg:4269 +proj=longlat +ellps=GRS80 ",
                  "+datum=NAD83 +no_defs +towgs84=0,0,0")
LATLONG <- sp::CRS(LATLONG)
ALBEA <- paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 ",
                "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
ALBEA <- sp::CRS(ALBEA)

DD <- data.frame(site_no=sitefile$site_no,
                 dec_lat_va=sitefile$dec_lat_va,
                 dec_long_va=sitefile$dec_long_va,
                 CDA=log10(CDA), stringsAsFactors=FALSE)
DD <- merge(FDC, DD, all=TRUE)


DD <- SpatialPointsDataFrame(cbind(DD$dec_long_va, DD$dec_lat_va), DD,
                             proj4string=LATLONG)
DD <- spTransform(DD, ALBEA)
XY <- coordinates(DD)
DD$east <- XY[,1]/1000; DD$north <- XY[,2]/1000; rm(XY)

SO <- over(DD, spDNI_1998to2009)
DD$ANN_DNI <- SO$ANN_DNI
DD$JAN <- SO$JAN; DD$FEB <- SO$FEB
DD$MAR <- SO$MAR; DD$APR <- SO$APR
DD$MAY <- SO$MAY; DD$JUN <- SO$JUN
DD$JUL <- SO$JUL; DD$AUG <- SO$AUG
DD$SEP <- SO$SEP; DD$OCT <- SO$OCT
DD$NOV <- SO$NOV; DD$DEC <- SO$DEC
rm(SO)

DD$x <- DD$east; DD$y <- DD$north
DDo <- DD

bnd <- list(x=bnd_poly_aea[,1]/1000, y=bnd_poly_aea[,2]/1000)
bnd <- list(bnd)
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
knots <- knots[-c(1, 2, 4, 8, 9, 12, 13, 21),]
x <- knots$x; y <- knots$y
knots <- data.frame(x=c(x, -200), y=c(y, 800)); rm(x,y)
points(knots$x, knots$y, pch=16, cex=1.1, col=4)

knots_pplo <- knots

length(DDo$site_no)
length(DD$site_no)

#dd <- slot(DD, name="data")
#for(i in c(4:5, 45:79))
#n     <- aggregate(DD$n,     by=list(DD$site_no), sum)$x
#nzero <- aggregate(DD$nzero, by=list(DD$site_no), sum)$x

# storages are in acre-ft
# 1 km2 = 247.104393047 acres
DD$flood_storage <- (DD$acc_nid_storage - DD$acc_norm_storage)/(DD$acc_basin_area*247.104393047)
DD[DD$flood_storage < 0,] # two sites: 02295420 and 02296750
DD$flood_storage <- abs(DD$flood_storage)
DD$flood_storage <- log10(DD$flood_storage+.01)
plot(qnorm(pp(DD$flood_storage)), sort(DD$flood_storage))
message("Maximum log10offets of flood_storage=",max(DD$flood_storage))


DD$ppt_mean        <- log10(DD$ppt_mean)
DD$temp_mean       <- log10(DD$temp_mean)
DD$acc_basin_area  <- log10(DD$acc_basin_area)
DD$acc_basin_slope <- log10(DD$acc_basin_slope)
plot(DD$CDA, DD$acc_basin_area, lwd=0.5,
     xlab="log10(NWIS CDA)", ylab="log10(NHDplus basin area")
abline(0,1)
abline(1/3,1,lty=2); abline(-1/3,1,lty=2)
abline(1/2,1,lty=2); abline(-1/2,1,lty=3)
mtext("Diagnostic check on watershed areas")
sites_of_area_bust <- unique(DD$site_no[abs(DD$acc_basin_area - DD$CDA) > 1/2])
for(site in sites_of_area_bust) {
  points(DD$CDA[DD$site_no == site], DD$acc_basin_area[DD$site_no == site], pch=16, col=2)
  DD <- DD[DD$site_no != site,]
}
text(0,5, paste(sites_of_area_bust, collapse=", "), cex=0.6, pos=4)
#points(DD$CDA[DD$site_no == "08167000"],
#       DD$acc_basin_area[DD$site_no == "08167000"], pch=16, col=4)

DD$decade   <- as.factor(DD$decade)
DD$bedperm  <- as.factor(DD$bedperm)
DD$aquifers <- as.factor(DD$aquifers)
DD$soller   <- as.factor(DD$soller)
DD$hlr      <- as.factor(DD$hlr)
DD$ecol3    <- as.factor(DD$ecol3)
DD$physio   <- as.factor(DD$physio)
DD$statsgo  <- as.factor(DD$statsgo)

DD$barren <- 2*asin(sqrt(DD$barren/100))
DD$cultivated_cropland <- 2*asin(sqrt(DD$cultivated_cropland/100))
DD$deciduous_forest    <- 2*asin(sqrt(DD$deciduous_forest/100))
DD$developed           <- 2*asin(sqrt(DD$developed/100))
DD$evergreen_forest    <- 2*asin(sqrt(DD$evergreen_forest/100))
DD$grassland           <- 2*asin(sqrt(DD$grassland/100))
DD$hay_pasture         <- 2*asin(sqrt(DD$hay_pasture/100))
DD$herbaceous_wetland  <- 2*asin(sqrt(DD$herbaceous_wetland/100))
DD$mixed_forest        <- 2*asin(sqrt(DD$mixed_forest/100))
# Note perennial_ice_snow is 0.00 throughout
DD$perennial_ice_snow  <- 2*asin(sqrt(DD$perennial_ice_snow/100))
DD$shrubland           <- 2*asin(sqrt(DD$shrubland/100))
DD$water               <- 2*asin(sqrt(DD$water/100))
DD$woody_wetland       <- 2*asin(sqrt(DD$woody_wetland/100))


DD$alt_ecol3 <- "-0"
#DD$alt_ecol3[DD$ecol3 == "ecol3_26"] <- "-26"
DD$alt_ecol3[DD$ecol3 == "ecol3_27"] <- "-27"
#DD$alt_ecol3[DD$ecol3 == "ecol3_30"] <- "-30"
DD$alt_ecol3[DD$ecol3 == "ecol3_31"] <- "-31"
#DD$alt_ecol3[DD$ecol3 == "ecol3_32"] <- "-32"
#DD$alt_ecol3[DD$ecol3 == "ecol3_33"] <- "-33"
DD$alt_ecol3[DD$ecol3 == "ecol3_34"] <- "-34"
#DD$alt_ecol3[DD$ecol3 == "ecol3_36"] <- "-36"
DD$alt_ecol3[DD$ecol3 == "ecol3_75"] <- "-75"
DD$alt_ecol3 <- as.factor(DD$alt_ecol3)

DD$trimmed_aquifers <- "aother"
DD$trimmed_aquifers[DD$aquifers == "TOT_AQ111"] <- "TOT_AQ111(surficial aquifer system)"
DD$trimmed_aquifers[DD$aquifers == "TOT_AQ413"] <- "TOT_AQ413(Floridan aquifer system)"
DD$trimmed_aquifers[DD$aquifers == "TOT_AQ201"] <- "TOT_AQ201(coastal lowlands aquifer system)"
DD$trimmed_aquifers <- as.factor(DD$trimmed_aquifers)

# This is the coastal plain
DD$physio <- relevel(DD$physio, "cat_physio_3")
DD$alt_physio <- "other"
DD$alt_physio[DD$physio == "cat_physio_3"] <- "Coastal Plain"
DD$alt_physio[DD$physio == "cat_physio_8"] <- "Appalachian Plateaus"
DD$alt_physio[DD$physio == "cat_physio_12"] <- "Central Lowland"
DD$alt_physio[DD$physio == "cat_physio_13"] <- "Great Plains"
DD$alt_physio <- as.factor(DD$alt_physio)
DD$alt_physio <- relevel(DD$alt_physio, "other")



DD$edwards_rechzone <- 0
DD$edwards_rechzone[DD$site_no == "08156800"] <- 1
DD$edwards_rechzone[DD$site_no == "08155300"] <- 1
DD$edwards_rechzone[DD$site_no == "08155400"] <- 1
DD$edwards_rechzone[DD$site_no == "08181400"] <- 1
DD$edwards_rechzone[DD$site_no == "08184000"] <- 1
DD$edwards_rechzone[DD$site_no == "08185000"] <- 1
DD$edwards_rechzone[DD$site_no == "08190500"] <- 1
DD$edwards_rechzone[DD$site_no == "08197500"] <- 1
DD$edwards_rechzone[DD$site_no == "08198500"] <- 1
DD$edwards_rechzone[DD$site_no == "08200700"] <- 1
DD$edwards_rechzone[DD$site_no == "08202700"] <- 1
DD$edwards_rechzone <- as.factor(DD$edwards_rechzone)

D <- DD;

D <- D[D$edwards_rechzone != "1",]

save(D, DD, DDo, knots, knots_pplo, bnd, file="DEMO.RData")

duan_smearing_estimator <- function(model) { sum(10^residuals(model))/length(residuals(model)) }


# [45] "ppt_mean"            "ppt_sd"              "temp_mean"           "temp_sd"
# [49] "tot_hdens"           "tot_major"           "tot_ndams"           "tot_nid_storage"
# [53] "tot_norm_storage"    "barren"              "cultivated_cropland" "deciduous_forest"
# [57] "developed"           "evergreen_forest"    "grassland"           "hay_pasture"
# [61] "herbaceous_wetland"  "mixed_forest"        "perennial_ice_snow"  "shrubland"
# [65] "water"               "woody_wetland"       "tot_bfi"             "sinuosity"
# [69] "length_km"           "area_sqkm"           "strm_dens"           "tot_twi"
# [73] "acc_basin_area"      "acc_basin_slope"     "tot_elev_mean"       "tot_elev_min"
# [77] "tot_elev_max"        "tot_total_road_dens" "tot_rdx"             "bedperm"
# [81] "aquifers"            "soller"              "hlr"                 "ecol3"
# [85] "physio"              "statsgo"

# It appears critical that the boundary have variables named say x,y
# The knots have the same names (x,y) and most difficult to figure out
# the x, y must be passed as same names in to the s(...). The x,y is not
# the critical piece, it is that they are all the same literal string.
# For example, v,w would work too.
plot(bnd[[1]]$x,bnd[[1]]$y,type="l", col=8, lwd=.6)
points(D$east[D$nzero == 0], D$north[D$nzero == 0], pch=4, lwd=.5, cex=0.9, col=rgb(0,.5,0.5,.5))
points(D$east[D$nzero > 0 & D$nzero <= 800], D$north[D$nzero > 0 & D$nzero <= 800], pch=4, lwd=.5, cex=0.9, col=rgb(1,0,0.5,.5))
points(D$east[D$nzero > 3653-2000], D$north[D$nzero > 3653-2000], pch=16, lwd=.5, cex=0.9, col=rgb(0.5,0,1))


family <- "gaussian"
Z <- D
x <- Z$east; y <- Z$north
Z$flowtime <- log10(Z$n - Z$nzero)
Zc <- Surv(Z$flowtime, Z$nzero != 0, type="right")
SM1 <- survreg(Zc~acc_basin_area+ppt_mean+temp_mean+acc_basin_slope+flood_storage+developed+ANN_DNI+bedperm+decade-1, data=Z, dist=family)
P1 <- P1o <- predict(SM1); P1[P1 > log10(3653)] <- log10(3653)
plot(Z$flowtime, P1, xlim=c(2.3,3.6), ylim=c(3,3.6), pch=16, col=rgb(0,0,1,.2))
abline(0,1)

SM0 <- survreg(Zc~acc_basin_area+
                    ppt_mean+temp_mean+ANN_DNI+
                    developed+grassland+
                    bedperm+decade-1, data=Z, dist=family)
Psurv0 <- Psurv0o <- predict(SM0); Psurv0[Psurv0 > log10(3653)] <- log10(3653)
plot(Z$flowtime, Psurv0, xlim=c(2.3,3.6), ylim=c(3,3.6), pch=16, col=rgb(0,0,1,.2))
abline(0,1)



library(cenGAM)
Z <- D
x <- Z$east; y <- Z$north
Z$left.threshold <- rep(0, length(Z$nzero))
Z$right.threshold <- log10(Z$n)
Z$flowtime <- log10(Z$n - Z$nzero)
GM0 <- gam(flowtime~acc_basin_area+
                    ppt_mean+s(temp_mean)+ANN_DNI+
                    developed+s(grassland)+
                    bedperm+decade-1,
           family=tobit1(left.threshold=  Z$left.threshold,
                         right.threshold=Z$right.threshold), data=Z); summary(GM0)
P0 <- P0o <- predict(GM0); P0[P0 > log10(3653)] <- log10(3653)
plot(Psurv0, P0, xlab="Survival Regression: log10(days of decadal streamflow)",
                 ylab="Censored GAM: log10(days of decadal streamflow)")
abline(0,1)
mtext("Model structure not quite exact because of convergence issues in GAM")
# WHA I would rather not have the smooth on temp_mean and grassland but these
# appear needed to get the GAM to converge. So the model does not quite match
# that from survreg().

GM1 <- gam(flowtime~acc_basin_area+
                    s(ppt_mean, k=5)+s(temp_mean)+s(ANN_DNI)+
                    developed+s(grassland)+
                    bedperm+decade-1,
           family=tobit1(left.threshold=  Z$left.threshold,
                         right.threshold=Z$right.threshold), data=Z); summary(GM1)
GM2 <- gam(flowtime~acc_basin_area+
                    s(ppt_mean, k=5)+s(temp_mean, k=4)+s(ANN_DNI, k=7)+
                    developed+s(grassland)+
                    bedperm+decade-1+
        s(x,y, bs="so", xt=list(bnd=bnd)), knots=knots,
           family=tobit1(left.threshold=  Z$left.threshold,
                         right.threshold=Z$right.threshold), data=Z); summary(GM2)

P <- Po <- predict(SM1); P[P > log10(3653)] <- log10(3653)
C1 <- C1o <- predict(GM1); C1[C1 > log10(3653)] <- log10(3653)
C2 <- C2o <- predict(GM2); C2[C2 > log10(3653)] <- log10(3653)
pa <- abs(P - Z$flowtime)
ca <- abs(C1- Z$flowtime)
cb <- abs(C2- Z$flowtime)
summary(pa)
summary(ca)
summary(cb)


plot(Z$flowtime, P, xlim=c(2.9,3.6), ylim=c(3,3.6), type="n")
abline(0,1)
points(Z$flowtime ,P, pch=16, col=rgb(0,0,1,.5), cex=0.5)
points(Z$flowtime,C2, pch=16, col=rgb(1,0,0,.2), cex=0.5)
#for(i in 1:length(P)) {
#  if(abs(P[i]-C1[i]) < .01) next
#  try(arrows(Z$flowtime[i],P[i],Z$flowtime[i],C1[i],
#             lwd=0.5, angle=10, length=.1, col=4), silent=TRUE)
#}
for(i in 1:length(P)) {
  if(abs(P[i]-C2[i]) < .02) next
  try(arrows(Z$flowtime[i],P[i],Z$flowtime[i],C2[i],
             lwd=0.5, angle=20, length=.1, col=2), silent=TRUE)
}

plot( qnorm(pp(pa)), sort(pa), type="l", xlim=c(0,4))
lines(qnorm(pp(ca)), sort(ca), col=4)
lines(qnorm(pp(cb)), sort(cb), col=2)
points(qnorm(pp(cb)), sort(cb), col=2, lwd=0.5)

Ppplo  <- (3653-10^P) /3653
C1pplo <- (3653-10^C1)/3653
C2pplo <- (3653-10^C2)/3653
plot(  Z$pplo, Ppplo, xlim=c(0,1), ylim=c(0,1))
points(Z$pplo, C2pplo, col=2)
points(Z$pplo, C2pplo, col=4)
abline(0,1)

PPLO <- GM2
save(D, SM0, SM1, GM1, GM2, PPLO, file="PPLOS.RData")



plot(bnd[[1]]$x,bnd[[1]]$y,type="l", col=8, lwd=.6)
points(D$east[D$nzero == 0], D$north[D$nzero == 0], pch=4, lwd=.5, cex=0.9, col=rgb(0,.5,0.5,.5))
points(D$east[D$nzero > 0 & D$nzero <= 800], D$north[D$nzero > 0 & D$nzero <= 800], pch=4, lwd=.5, cex=0.9, col=rgb(1,0,0.5,.5))
points(D$east[D$nzero > 3653-2000], D$north[D$nzero > 3653-2000], pch=16, lwd=.5, cex=0.9, col=rgb(0.5,0,1))
mtext("Raw of PPLO")



plot(bnd[[1]]$x,bnd[[1]]$y,type="l", col=8, lwd=.6)
points(D$east[P == log10(3653)], D$north[P == log10(3653)], pch=4, lwd=.5, cex=0.9, col=rgb(0,.5,0.5,.5))
points(D$east[P < log10(3653)], D$north[P < log10(3653)], pch=4, lwd=.5, cex=0.9, col=rgb(1,0,0.5,.5))
points(D$east[P < log10(2000)], D$north[P < log10(2000)], pch=16, lwd=.5, cex=0.9, col=rgb(0.5,0,1,.5))
mtext("Predictions of PPLO")


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

###### END PPLO

# [45] "ppt_mean"            "ppt_sd"              "temp_mean"           "temp_sd"
# [49] "tot_hdens"           "tot_major"           "tot_ndams"           "tot_nid_storage"
# [53] "tot_norm_storage"    "barren"              "cultivated_cropland" "deciduous_forest"
# [57] "developed"           "evergreen_forest"    "grassland"           "hay_pasture"
# [61] "herbaceous_wetland"  "mixed_forest"        "perennial_ice_snow"  "shrubland"
# [65] "water"               "woody_wetland"       "tot_bfi"             "sinuosity"
# [69] "length_km"           "area_sqkm"           "strm_dens"           "tot_twi"
# [73] "acc_basin_area"      "acc_basin_slope"     "tot_elev_mean"       "tot_elev_min"
# [77] "tot_elev_max"        "tot_total_road_dens" "tot_rdx"             "bedperm"
# [81] "aquifers"            "soller"              "hlr"                 "ecol3"
# [85] "physio"              "statsgo"

Z <- D
x <- Z$x; y <- Z$y
z <- log10(Z$L1)
L1   <- gam(z~acc_basin_area +
              s(ppt_mean, k=5) + s(temp_mean, k=4) + s(ANN_DNI, k=7)+
              developed+
              mixed_forest+shrubland+
              decade-1+
              s(x, y, bs="so", xt=list(bnd=bnd)), knots=knots, data=Z,
              family="gaussian")
L1$duan_smearing <- duan_smearing_estimator(L1)
pdf("L1.pdf", useDingbats=FALSE)
  plot(z, fitted.values(L1))
  abline(0,1)
  plot(L1, scheme=2)
  points(Z$x, Z$y, pch=4, lwd=.5, cex=0.9, col=8)
  points(knots$x, knots$y, pch=16, cex=1.1, col=4)
  text(100, 500, "L1")
dev.off()

z <- D$L2/D$L1 # --------------------------- Coefficient of L-variation
T2   <- gam(z~acc_basin_area +
              s(ppt_mean, k=5) + s(temp_mean, k=4) + s(ANN_DNI, k=7)+
              developed+
              mixed_forest+shrubland+
              decade-1+
              s(x, y, bs="so", xt=list(bnd=bnd)), knots=knots, data=D,
              family="gaussian")
pdf("T2.pdf", useDingbats=FALSE)
  plot(z, fitted.values(T2))
  abline(0,1)
  plot(T2, scheme=2)
  points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
  points(knots$x, knots$y, pch=16, cex=1.1, col=4)
  text(100, 500, "Tau2")
dev.off()

z <- D$T3      # --------------------------- L-skew
T3   <- gam(z~acc_basin_area +
              s(ppt_mean, k=5) + s(temp_mean, k=4) + s(ANN_DNI, k=7)+
              developed+
              mixed_forest+shrubland+
              decade-1+
              s(x, y, bs="so", xt=list(bnd=bnd)), knots=knots, data=D,
              family="gaussian")
pdf("T3.pdf", useDingbats=FALSE)
  plot(z, fitted.values(T3))
  abline(0,1)
  plot(T3, scheme=2)
  points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
  points(knots$x, knots$y, pch=16, cex=1.1, col=4)
  text(100, 500, "Tau3")
dev.off()


z <- D$T4      # --------------------------- L-kurtosis
T4   <- gam(z~acc_basin_area +
              s(ppt_mean, k=5) + s(temp_mean, k=4) + s(ANN_DNI, k=7)+
              developed+
              mixed_forest+shrubland+
              decade-1+
              s(x, y, bs="so", xt=list(bnd=bnd)), knots=knots, data=D,
              family="gaussian")
pdf("T4.pdf", useDingbats=FALSE)
  plot(z, fitted.values(T4))
  abline(0,1)
  plot(T4, scheme=2)
  points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
  points(knots$x, knots$y, pch=16, cex=1.1, col=4)
  text(100, 500, "Tau4")
dev.off()



z <- D$T5      # --------------------------- Fifth L-moment ratio
T5   <- gam(z~acc_basin_area +
              s(ppt_mean, k=5) + s(temp_mean, k=4) + s(ANN_DNI, k=7)+
              developed+
              mixed_forest+shrubland+
              decade-1+
              s(x, y, bs="so", xt=list(bnd=bnd)), knots=knots, data=D,
              family="gaussian")
pdf("T5.pdf", useDingbats=FALSE)
  plot(z, fitted.values(T5))
  abline(0,1)
  plot(T5, scheme=2)
  points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
  points(knots$x, knots$y, pch=16, cex=1.1, col=4)
  text(100, 500, "Tau5")
dev.off()


z <- D$T6      # --------------------------- Sixth L-moment ratio
T6   <- gam(z~acc_basin_area +
              s(ppt_mean, k=5) + s(temp_mean, k=4) + s(ANN_DNI, k=7)+
              developed+
              mixed_forest+shrubland+
              decade-1+
              s(x, y, bs="so", xt=list(bnd=bnd)), knots=knots, data=D,
              family="gaussian")
pdf("T6.pdf", useDingbats=FALSE)
  plot(z, fitted.values(T6))
  abline(0,1)
  plot(T6, scheme=2)
  points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
  points(knots$x, knots$y, pch=16, cex=1.1, col=4)
  text(100, 500, "Tau6")
dev.off()


z <- log10(D$f50+1)      # --------------------------- Sixth L-moment ratio
Q50   <- gam(z~acc_basin_area +
               s(ppt_mean, k=5) + s(temp_mean, k=4) + s(ANN_DNI, k=7)+
               developed+
               mixed_forest+shrubland+
               decade-1+
               s(x, y, bs="so", xt=list(bnd=bnd)), knots=knots, data=D,
             family="gaussian")
pdf("Q50.pdf", useDingbats=FALSE)
plot(z, fitted.values(Q50))
abline(0,1)
plot(Q50, scheme=2)
points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
points(knots$x, knots$y, pch=16, cex=1.1, col=4)
text(100, 500, "Q50")
dev.off()


sink("right_tail_flowing_fdc.txt")
z <- log10(D$f90+1)      # --------------------------- Sixth L-moment ratio
Q90   <- gam(z~acc_basin_area  + s(flood_storage, k=7) +
              s(ppt_mean, k=5) + s(ANN_DNI, k=7)+
              developed+
              mixed_forest+shrubland+
              decade-1+
              s(x, y, bs="so", xt=list(bnd=bnd)), knots=knots, data=D,
              family="gaussian")
print(summary(Q90))
pdf("Q90.pdf", useDingbats=FALSE)
  plot(z, fitted.values(Q90))
  abline(0,1)
  plot(Q90, scheme=2, residuals=TRUE)
  points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
  points(knots$x, knots$y, pch=16, cex=1.1, col=4)
  text(100, 500, "Q90")
dev.off()


z <- log10(D$f95+1)      # --------------------------- Sixth L-moment ratio
Q95   <- gam(z~acc_basin_area + s(flood_storage, k=7) +
               s(ppt_mean, k=5) + s(ANN_DNI, k=7)+
               developed+
               mixed_forest+shrubland+
               decade-1+
               s(x, y, bs="so", xt=list(bnd=bnd)), knots=knots, data=D,
             family="gaussian")
print(summary(Q95))
pdf("Q95.pdf", useDingbats=FALSE)
plot(z, fitted.values(Q95))
abline(0,1)
plot(Q95, scheme=2, residuals=TRUE)
points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
points(knots$x, knots$y, pch=16, cex=1.1, col=4)
text(100, 500, "Q95")
dev.off()


z <- log10(D$f98+1)      # --------------------------- Sixth L-moment ratio
Q98   <- gam(z~acc_basin_area + s(flood_storage, k=7) +
               s(ppt_mean, k=5) + s(ANN_DNI, k=7)+
               developed+
               mixed_forest+shrubland+
               decade-1+
               s(x, y, bs="so", xt=list(bnd=bnd)), knots=knots, data=D,
             family="gaussian")
print(summary(Q98))
pdf("Q98.pdf", useDingbats=FALSE)
plot(z, fitted.values(Q98))
abline(0,1)
plot(Q98, scheme=2, residuals=TRUE)
points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
points(knots$x, knots$y, pch=16, cex=1.1, col=4)
text(100, 500, "Q98")
dev.off()

z <- log10(D$f99+1)      # --------------------------- Sixth L-moment ratio
Q99   <- gam(z~acc_basin_area + s(flood_storage, k=7) +
               s(ppt_mean, k=5) + s(ANN_DNI, k=7)+
               developed+
               mixed_forest+shrubland+
               decade-1+
               s(x, y, bs="so", xt=list(bnd=bnd)), knots=knots, data=D,
             family="gaussian")
print(summary(Q99))
pdf("Q99.pdf", useDingbats=FALSE)
plot(z, fitted.values(Q99))
abline(0,1)
plot(Q99, scheme=2, residuals=TRUE)
points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
points(knots$x, knots$y, pch=16, cex=1.1, col=4)
text(100, 500, "Q99")
dev.off()


z <- log10(D$f99.9+1)      # --------------------------- Sixth L-moment ratio
Q99p9   <- gam(z~acc_basin_area + s(flood_storage, k=7) +
               s(ppt_mean, k=5) + s(ANN_DNI, k=7)+
               developed+
               mixed_forest+shrubland+
               decade-1+
               s(x, y, bs="so", xt=list(bnd=bnd)), knots=knots, data=D,
             family="gaussian")
print(summary(Q99p9))
pdf("Q99p9.pdf", useDingbats=FALSE)
plot(z, fitted.values(Q99p9))
abline(0,1)
plot(Q99p9, scheme=2, residuals=TRUE)
points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
points(knots$x, knots$y, pch=16, cex=1.1, col=4)
text(100, 500, "Q99p9")
dev.off()




save(D, PPLO, L1, T2, T3, T4, T5, T6, Q50, Q90, Q95, Q98, Q99, Q99p9, file="Models.RData")
