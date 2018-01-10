library(feather)
library(proxy)
library(akqdecay)
library(mgcv)
library(sp)
library(rgeos)
library(lmomco)
library(Lmoments)

FDC <- read_feather(file.choose()) # "all_gage_data.feather"
load(file.choose()) # "spRESTORE_MGCV_BND.RData"
load(file.choose()) # spDNI_1998to2009.RData

sites <- unique(FDC$site_no)
sitefile <- dataRetrieval::readNWISsite(sites)
sitefile <- sitefile[sitefile$agency_cd != "USCE",]

CDA <- sitefile$contrib_drain_area_va
CDA[is.na(CDA)] <- sitefile$drain_area_va[is.na(CDA)]
sitefile$contrib_drain_area_va <- CDA


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
DD <- DD[complete.cases(DD),]

DD <- SpatialPointsDataFrame(cbind(DD$dec_long_va, DD$dec_lat_va), DD,
                            proj4string=LATLONG)
DD <- spTransform(DD, ALBEA)
XY <- coordinates(DD)
DD$east <- XY[,1]/1000; DD$north <- XY[,2]/1000; rm(XY)

DD$prob_has_noflow <- 0
DD$prob_has_noflow[DD$nzero > 0] <- 1


SO <- over(DD, spDNI_1998to2009)
DD$ANN_DNI <- SO$ANN_DNI
DD$JAN <- SO$JAN; DD$FEB <- SO$FEB
DD$MAR <- SO$MAR; DD$APR <- SO$APR
DD$MAY <- SO$MAY; DD$JUN <- SO$JUN
DD$JUL <- SO$JUL; DD$AUG <- SO$AUG
DD$SEP <- SO$SEP; DD$OCT <- SO$OCT
DD$NOV <- SO$NOV; DD$DEC <- SO$DEC
rm(SO)

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
knots <- knots[-c(3, 13, 10, 11, 12, 29, 21, 22, 27, 16, 6, 5),]
x <- knots$x; y <- knots$y
knots <- data.frame(x=c(x, 1375), y=c(y, 600)); rm(x,y)
points(knots$x, knots$y, pch=16, cex=1.1, col=4)

knots_pplo <- knots[-17, ]

# The as.numeric() is needed head of pplo use later
plogit <- function(eta) as.numeric(exp(eta)/(exp(eta)+1))
duan_smearing_estimator <- function(model) { sum(10^residuals(L1))/length(residuals(L1)) }

DD$x <- DD$east; DD$y <- DD$north
DD$flood_storage <- (DD$tot_nid_storage - DD$tot_norm_storage)/10^DD$CDA
DD[DD$flood_storage < 0,] # two sites: 02295420 and 02296750
DD$flood_storage[DD$flood_storage < 0] <- (DD$tot_norm_storage[DD$flood_storage < 0] -
                                            DD$tot_nid_storage[DD$flood_storage < 0]) /
                                                     10^DD$CDA[DD$flood_storage < 0]



DD$tot_basin_slope <- log10(DD$tot_basin_slope)
DD$decade <- as.factor(DD$decade)
DD$aquifers <- as.factor(DD$aquifers)

D <- DD;

# [45] "ppt_mean"            "ppt_sd"              "temp_mean"           "temp_sd"
# [49] "tot_hdens"           "tot_major"           "tot_ndams"           "tot_nid_storage"
# [53] "tot_norm_storage"    "barren"              "cultivated_cropland" "deciduous_forest"
# [57] "developed"           "evergreen_forest"    "grassland"           "hay_pasture"
# [61] "herbaceous_wetland"  "mixed_forest"        "perennial_ice_snow"  "shrubland"
# [65] "water"               "woody_wetland"       "tot_bfi"             "sinuosity"
# [69] "length_km"           "area_sqkm"           "strm_dens"           "tot_twi"
# [73] "tot_basin_area"      "tot_basin_slope"     "tot_elev_mean"       "tot_elev_min"
# [77] "tot_elev_max"        "tot_total_road_dens" "tot_rdx"             "bedperm"
# [81] "aquifers"            "soller"              "hlr"                 "ecol3"
# [85] "physio"              "statsgo"

# It appears critical that the boundary have variables named say x,y
# The knots have the same names (x,y) and most difficult to figure out
# the x, y must be passed as same names in to the s(...). The x,y is not
# the critical piece, it is that they are all the same literal string.
# For example, v,w would work too.
Z <- D[D$nzero > 0,]
#x <- Z$east; y <- Z$north # Please see the extraction and setting to x,y
#CDA <- Z$CDA; ANN_DNI <- Z$ANN_DNI; MAY <- Z$MAY; DEC <- Z$DEC
Z$z <- D$pplo; Z$mf <- mf
PPLO <- gam(z~s(CDA, bs="cr", k=6)+s(ANN_DNI, bs="cr", k=6)+
              te(ppt_mean, temp_mean)+mixed_forest+developed+
              s(MAY, DEC, bs="tp", k=6)+
              s(x, y, bs="so", xt=list(bnd=bnd)),
              knots=knots_pplo, data=Z, family="quasibinomial")
pdf("PPLO.pdf", useDingbats=FALSE)
  plot(Z$ANN_DNI, Z$z, pch=16, col=8)
  points(Z$ANN_DNI, fitted.values(PPLO))
  plot(PPLO, scheme=2)
  points(Z$x, Z$y, pch=4, lwd=.5, cex=0.9, col=8)
  points(knots_pplo$x, knots_pplo$y, pch=16, cex=1.1, col=4)
  text(100, 500, "Probability level of noflows", pos=4)
dev.off()


Z <- D
z <- Z$prob_has_noflow
ZERO <- gam(z~s(CDA, bs="cr", k=6)+s(ANN_DNI, bs="cr", k=6)+
              te(ppt_mean, temp_mean)+
              s(MAY, DEC, bs="tp", k=6)+
              s(x, y, bs="so", xt=list(bnd=bnd))+
              developed,
              knots=knots_pplo, data=Z, family="quasibinomial")
pdf("ZERO.pdf", useDingbats=FALSE)
  plot(Z$ANN_DNI, z, pch=16, col=8);
  points(Z$ANN_DNI, fitted.values(ZERO))
  plot(ZERO, scheme=2)
  points(x, y, pch=4, lwd=.5, cex=0.9, col=8)
  points(knots_pplo$x, knots_pplo$y, pch=16, cex=1.1, col=4)
  text(100, 500, "Probability model of site having at least one noflow", pos=4)
dev.off()


z <- log10(D$L1)
#L1   <- gam(z~s(CDA, bs="cr")+s(ANN_DNI, bs="cr")+
#              s(MAY, DEC, bs="tp", k=5)+
#              s(x, y, bs="so", xt=list(bnd=bnd)), knots=knots,
#              family="gaussian")
#L1   <- gam(z~CDA*ANN_DNI+
#              s(MAY, DEC, bs="tp", k=5)+
#              s(x, y, bs="so", xt=list(bnd=bnd)), knots=knots,
#              family="gaussian")
L1   <- gam(z~te(CDA, ANN_DNI, bs="tp")+
              s(MAY, DEC, bs="tp", k=5)+
              s(x, y, bs="so", xt=list(bnd=bnd)), knots=knots,
              family="gaussian")
L1$duan_smearing <- duan_smearing_estimator(L1)
pdf("L1.pdf", useDingbats=FALSE)
  plot(z, fitted.values(L1))
  abline(0,1)
  plot(L1, scheme=2)
  points(x, y, pch=4, lwd=.5, cex=0.9, col=8)
  points(knots$x, knots$y, pch=16, cex=1.1, col=4)
  text(100, 500, "L1")
dev.off()

z <- D$L2/D$L1 # --------------------------- Coefficient of L-variation
T2   <- gam(z~s(CDA, bs="cr")+s(ANN_DNI, bs="cr")+
              s(MAY, DEC, bs="tp", k=6)+
              s(x, y, bs="so", xt=list(bnd=bnd)), knots=knots,
              family="gaussian")
pdf("T2.pdf", useDingbats=FALSE)
  plot(z, fitted.values(T2))
  abline(0,1)
  plot(T2, scheme=2)
  points(x, y, pch=4, lwd=.5, cex=0.9, col=8)
  points(knots$x, knots$y, pch=16, cex=1.1, col=4)
  text(100, 500, "Tau2")
dev.off()

z <- D$T3      # --------------------------- L-skew
T3   <- gam(z~s(CDA, bs="cr")+s(ANN_DNI, bs="cr")+
              s(MAY, DEC, bs="tp", k=6)+
              s(x, y, bs="so", xt=list(bnd=bnd)), knots=knots,
              family="gaussian")
pdf("T3.pdf", useDingbats=FALSE)
  plot(z, fitted.values(T3))
  abline(0,1)
  plot(T3, scheme=2)
  points(x, y, pch=4, lwd=.5, cex=0.9, col=8)
  points(knots$x, knots$y, pch=16, cex=1.1, col=4)
  text(100, 500, "Tau3")
dev.off()


z <- D$T4      # --------------------------- L-kurtosis
T4   <- gam(z~s(CDA, bs="cr")+s(ANN_DNI, bs="cr")+
              s(MAY, DEC, bs="tp", k=6)+
              s(x, y, bs="so", xt=list(bnd=bnd)), knots=knots,
              family="gaussian")
pdf("T4.pdf", useDingbats=FALSE)
  plot(z, fitted.values(T4))
  abline(0,1)
  plot(T4, scheme=2)
  points(x, y, pch=4, lwd=.5, cex=0.9, col=8)
  points(knots$x, knots$y, pch=16, cex=1.1, col=4)
  text(100, 500, "Tau4")
dev.off()



z <- D$T5      # --------------------------- Fifth L-moment ratio
T5   <- gam(z~s(CDA, bs="cr")+s(ANN_DNI, bs="cr")+
              s(MAY, DEC, bs="tp", k=6)+
              s(x, y, bs="so", xt=list(bnd=bnd)), knots=knots,
              family="gaussian")
pdf("T5.pdf", useDingbats=FALSE)
  plot(z, fitted.values(T5))
  abline(0,1)
  plot(T5, scheme=2)
  points(x, y, pch=4, lwd=.5, cex=0.9, col=8)
  points(knots$x, knots$y, pch=16, cex=1.1, col=4)
  text(100, 500, "Tau5")
dev.off()


z <- D$T6      # --------------------------- Sixth L-moment ratio
T6   <- gam(z~s(CDA, bs="cr")+s(ANN_DNI, bs="cr")+
              s(MAY, DEC, bs="tp", k=6)+
              s(x, y, bs="so", xt=list(bnd=bnd)), knots=knots,
              family="gaussian")
pdf("T6.pdf", useDingbats=FALSE)
  plot(z, fitted.values(T6))
  abline(0,1)
  plot(T6, scheme=2)
  points(x, y, pch=4, lwd=.5, cex=0.9, col=8)
  points(knots$x, knots$y, pch=16, cex=1.1, col=4)
  text(100, 500, "Tau6")
dev.off()


#--***---***---***---***---***---***---***---***---***---***---***---***---***
#                     DEMONSTRATION OF UNGAGED ESTIMATION
#--***---***---***---***---***---***---***---***---***---***---***---***---***
newdata <- data.frame(x=-50, y=1000, CDA=2.9, ANN_DNI=4.9, MAY=4.8, DEC=4.1)

prob_has_noflow <- plogit(predict(ZERO, newdata=newdata))
pplo <- plogit(predict(PPLO, newdata=newdata))
if(prob_has_noflow < 0.5) pplo <- 0

l1 <- predict(L1, newdata=newdata); l1 <- 10^l1
t2 <- predict(T2, newdata=newdata)
t3 <- predict(T3, newdata=newdata)
t4 <- predict(T4, newdata=newdata)
t5 <- predict(T5, newdata=newdata)
t6 <- predict(T6, newdata=newdata)
lmrvec <- c(l1, l1*t2, t3, t4, t5, t6)
lmr <- vec2lmom(lmrvec)
par <- lmom2par(lmr, type="pe3")

db <- data.frame(x=x, y=y, CDA=CDA, ANN_DNI=ANN_DNI, MAY=MAY, DEC=DEC)
dists <- proxy::dist(db, newdata)
site <- D$site[(1:length(dists))[dists == min(dists)]]
fdc <- D[D$site == site,]
fdc <- as.data.frame(fdc[,7:33])
fdc <- t(fdc); n <- length(fdc); fdc <- fdc[-c(n-1,n),]
ffs <- names(fdc); ffs <- as.numeric(gsub("f", "", ffs))/100

FF <- c(0.0001,seq(0.001,.999,by=.001),0.9999) # convenient nonexceed prob
nFF <- lmomco::f2flo(FF, pp=pplo)
qFF <- qnorm(lmomco::f2f(FF, pp=pplo))
ylim <- range(c(fdc, qlmomco(nFF, par)))
ylim[ylim[1] <= 0] <- c(0.01, ylim[2])
xlim <- qnorm(range(FF))

qdf <- qlmomco(nFF, par)

plot(qFF, qdf, type="l", log="y", xaxt="n",
     xlab="", ylab="QUANTILE, CFS", xlim=xlim, ylim=ylim, col=2, lwd=1.0)
lines(qFF, qdf*L1$duan_smearing, col=4, lwd=1.2)
add.lmomco.axis(las=2, tcl=0.5, side.type="NPP", npp.as.aep=TRUE) # x-axis
mtext("Ungaged FDC estimation (eight regionalization models involved)")
points(qnorm(ffs), fdc)

legend(-4,10000, c("Biased generalized normal",
                   "Unbiased (Duan smearing) generalized normal",
                   "WHA's moonshot (six L-moment fit)"), bty="n", lwd=c(1,1.2,1.4), col=c(2,4,3))

