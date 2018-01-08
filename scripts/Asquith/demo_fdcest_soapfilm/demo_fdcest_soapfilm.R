library(feather)
pplo <- read_feather(file.choose()) # "fdc_lmr_pplo.feather" or "fdc_lmr_pplo2010-15.feather"
# I am using fdc_lmr_pplo2010-15.feather for the testing here.

library(akqdecay)
library(mgcv)
library(sp)
library(rgeos)

sites <- unique(pplo$site)
sitefile <- dataRetrieval::readNWISsite(sites)
sitefile <- sitefile[sitefile$agency_cd != "USCE",]

CDA <- sitefile$contrib_drain_area_va
CDA[is.na(CDA)] <- sitefile$drain_area_va[is.na(CDA)]
sitefile$contrib_drain_area_va <- CDA

load(file.choose()) # "spRESTORE_MGCV_BND.RData"
LATLONG <- paste0("+init=epsg:4269 +proj=longlat +ellps=GRS80 ",
                  "+datum=NAD83 +no_defs +towgs84=0,0,0")
LATLONG <- sp::CRS(LATLONG)
ALBEA <- paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 ",
                "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
ALBEA <- sp::CRS(ALBEA)

D <- data.frame(site=sitefile$site_no,
                dec_lat_va=sitefile$dec_lat_va,
                dec_long_va=sitefile$dec_long_va,
                CDA=log10(CDA), stringsAsFactors=FALSE)
D <- merge(pplo, D, all=TRUE)
D <- D[complete.cases(D),]

D <- SpatialPointsDataFrame(cbind(D$dec_long_va,
                                  D$dec_lat_va), D,
                            proj4string=LATLONG)
D <- spTransform(D, ALBEA)
XY <- coordinates(D)
D$east <- XY[,1]/1000; D$north <- XY[,2]/1000

D$event <- 0
D$event[D$nzero > 0] <- 1

load(file.choose()) # spDNI_1998to2009.RData
SO <- over(D, spDNI_1998to2009)

D$ANN_DNI <- SO$ANN_DNI
D$JAN <- SO$JAN
D$FEB <- SO$FEB
D$MAR <- SO$MAR
D$APR <- SO$APR
D$MAY <- SO$MAY
D$JUN <- SO$JUN
D$JUL <- SO$JUL
D$AUG <- SO$AUG
D$SEP <- SO$SEP
D$OCT <- SO$OCT
D$NOV <- SO$NOV
D$DEC <- SO$DEC


fsb <- list(x=bnd_poly_aea[,1]/1000, y=bnd_poly_aea[,2]/1000)
fsb <- list(fsb)
plot(fsb[[1]]$x,fsb[[1]]$y,type="l", col=8, lwd=.6)
points(D$east, D$north, pch=3, lwd=.5, col=8)

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
ind <- mgcv::inSide(fsb,x,y)
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
points(knots$x, knots$y, pch=3, cex=0.6, col=4)

#knots <- data.frame(x=c(-250, 250, 500, 1000), y=rep(1000, 4))

# It appears critical that the boundary have variables named say x,y
# The knots have the same names (x,y) and most difficult to figure out
# the x, y must be passed as same names in to the s(...). The x,y is not
# the critical piece, it is that they are all the same literal string.
# For example, v,w would work too.
x <- D$east; y <- D$north # Please see the extraction and setting to x,y
CDA <- D$CDA; ANN_DNI <- D$ANN_DNI; MAY <- D$MAY; DEC <- D$DEC; z <- D$event
ZERO <- gam(z~s(CDA, bs="cr")+s(ANN_DNI, bs="cr")+
              s(x, y, bs="so", xt=list(bnd=fsb)), knots=knots,
              family="quasibinomial")
plot(ANN_DNI, z, pch=16, col=8);
points(ANN_DNI, fitted.values(ZERO))


Z <- D[D$nzero > 0,]
x <- Z$east; y <- Z$north
CDA <- Z$CDA; ANN_DNI <- Z$ANN_DNI; MAY <- Z$MAY; DEC <- Z$DEC; z <- Z$pplo
PPLO <- gam(z~s(CDA, bs="cr")+s(ANN_DNI, bs="cr")+
              s(MAY, DEC, bs="tp", k=6)+
              s(x, y, bs="so", xt=list(bnd=fsb)), knots=knots,
              family="quasibinomial")
plot(ANN_DNI, z, pch=16, col=8);
points(ANN_DNI, fitted.values(PPLO))



x <- D$east; y <- D$north # Please see the extraction and setting to x,y
CDA <- D$CDA; ANN_DNI <- D$ANN_DNI; MAY <- D$MAY; DEC <- D$DEC; z <- log10(D$L1)
#L1   <- gam(z~s(CDA, bs="cr")+s(ANN_DNI, bs="cr")+
#              s(MAY, DEC, bs="tp", k=5)+
#              s(x, y, bs="so", xt=list(bnd=fsb)), knots=knots,
#              family="gaussian")
#L1   <- gam(z~CDA*ANN_DNI+
#              s(MAY, DEC, bs="tp", k=5)+
#              s(x, y, bs="so", xt=list(bnd=fsb)), knots=knots,
#              family="gaussian")
L1   <- gam(z~te(CDA, ANN_DNI, bs="tp")+
              s(MAY, DEC, bs="tp", k=5)+
              s(x, y, bs="so", xt=list(bnd=fsb)), knots=knots,
              family="gaussian")
plot(z, fitted.values(L1))
abline(0,1)
duan_smearing_estimator <- function(model) { sum(10^residuals(L1))/length(residuals(L1)) }
duan_smearing <- duan_smearing_estimator(L1)


x <- D$east; y <- D$north # Please see the extraction and setting to x,y
CDA <- D$CDA; ANN_DNI <- D$ANN_DNI; MAY <- D$MAY; DEC <- D$DEC; z <- D$L2/D$L1
T2   <- gam(z~s(CDA, bs="cr")+s(ANN_DNI, bs="cr")+
              s(MAY, DEC, bs="tp", k=6)+
              s(x, y, bs="so", xt=list(bnd=fsb)), knots=knots,
              family="gaussian")
plot(z, fitted.values(T2))
abline(0,1)

x <- D$east; y <- D$north # Please see the extraction and setting to x,y
CDA <- D$CDA; ANN_DNI <- D$ANN_DNI; MAY <- D$MAY; DEC <- D$DEC; z <- D$T3
T3   <- gam(z~s(CDA, bs="cr")+s(ANN_DNI, bs="cr")+
              s(MAY, DEC, bs="tp", k=6)+
              s(x, y, bs="so", xt=list(bnd=fsb)), knots=knots,
              family="gaussian")
plot(z, fitted.values(T3))
abline(0,1)


x <- D$east; y <- D$north # Please see the extraction and setting to x,y
CDA <- D$CDA; ANN_DNI <- D$ANN_DNI; MAY <- D$MAY; DEC <- D$DEC; z <- D$T4
T4   <- gam(z~s(CDA, bs="cr")+s(ANN_DNI, bs="cr")+
              s(MAY, DEC, bs="tp", k=6)+
              s(x, y, bs="so", xt=list(bnd=fsb)), knots=knots,
              family="gaussian")
plot(z, fitted.values(T4))
abline(0,1)


x <- D$east; y <- D$north # Please see the extraction and setting to x,y
CDA <- D$CDA; ANN_DNI <- D$ANN_DNI; MAY <- D$MAY; DEC <- D$DEC; z <- D$T5
T5   <- gam(z~s(CDA, bs="cr")+s(ANN_DNI, bs="cr")+
              s(MAY, DEC, bs="tp", k=6)+
              s(x, y, bs="so", xt=list(bnd=fsb)), knots=knots,
              family="gaussian")
plot(z, fitted.values(T5))
abline(0,1)


x <- D$east; y <- D$north # Please see the extraction and setting to x,y
CDA <- D$CDA; ANN_DNI <- D$ANN_DNI; MAY <- D$MAY; DEC <- D$DEC; z <- D$T6
T6   <- gam(z~s(CDA, bs="cr")+s(ANN_DNI, bs="cr")+
              s(MAY, DEC, bs="tp", k=6)+
              s(x, y, bs="so", xt=list(bnd=fsb)), knots=knots,
              family="gaussian")
plot(z, fitted.values(T6))
abline(0,1)




##############################################
#  DEMONSTRATION OF UNGAGED ESTIMATION
##############################################
newdata <- data.frame(x=-50, y=1000, CDA=2.9, ANN_DNI=4.9, MAY=4.8, DEC=4.1)
eta <- predict(ZERO, newdata=newdata)
prob_has_noflow <- exp(eta)/(exp(eta)+1)

eta <- predict(PPLO, newdata=newdata)
pplo <- as.numeric(exp(eta)/(exp(eta)+1))
if(prob_has_noflow < 0.5) pplo <- 0

eta <- predict(L1, newdata=newdata)
l1 <- 10^eta
t2 <- predict(T2, newdata=newdata)
t3 <- predict(T3, newdata=newdata)
t4 <- predict(T4, newdata=newdata)
t5 <- predict(T5, newdata=newdata)
t6 <- predict(T6, newdata=newdata)
lmrvec <- c(l1, l1*t2, t3, t4, t5, t6)
lmr <- vec2lmom(lmrvec)
par <- lmom2par(lmr, type="gno")
FF <- c(0.0001,seq(0.001,.999,by=.001),0.9999) # convenience nonexceedance probabilities

db <- data.frame(x=x, y=y, CDA=CDA, ANN_DNI=ANN_DNI, MAY=MAY, DEC=DEC)
dists <- proxy::dist(db, newdata)
site <- D$site[(1:length(dists))[dists == min(dists)]]
fdc <- D[D$site == site,]
fdc <- as.data.frame(fdc[,7:33])
fdc <- t(fdc); n <- length(fdc); fdc <- fdc[-c(n-1,n),]
ffs <- names(fdc); ffs <- as.numeric(gsub("f", "", ffs))/100

ylim <- range(c(tmp, qlmomco(lmomco::f2flo(FF, pp=pplo), par)))
ylim[ylim[1] <= 0] <- c(0.01, ylim[2])
plot(qnorm(FF), lmomco::qlmomco(FF, par), type="l", log="y", xaxt="n",
     xlab="", ylab="QUANTILE, CFS", ylim=ylim, lty=2, lwd=0.8)
lines(qnorm(  lmomco::f2f(  FF, pp=pplo)),
      qlmomco(lmomco::f2flo(FF, pp=pplo), par), col=2, lwd=1)
lines(qnorm(  lmomco::f2f(  FF, pp=pplo)),
      qlmomco(lmomco::f2flo(FF, pp=pplo), par)*duan_smearing, col=4, lwd=1.2)
add.lmomco.axis(las=2, tcl=0.5, side.type="NPP", npp.as.aep=TRUE) # x-axis
mtext("Ungaged FDC estimation (dashed, naive; red, biased; green, bias corrected)")
points(qnorm(ffs), fdc)



