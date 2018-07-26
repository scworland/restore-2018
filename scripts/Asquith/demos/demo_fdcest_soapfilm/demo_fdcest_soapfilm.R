library(feather)
library(proxy)
library(akqdecay)
library(mgcv)
library(sp)
library(rgeos)

FDC <- read_feather(file.choose()) # "fdc_lmr_pplo.feather" or "fdc_lmr_pplo2010-15.feather"
# I am using fdc_lmr_pplo2010-15.feather for the testing here.
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

D <- data.frame(site_no=sitefile$site_no,
                dec_lat_va=sitefile$dec_lat_va,
                dec_long_va=sitefile$dec_long_va,
                CDA=log10(CDA), stringsAsFactors=FALSE)
D <- merge(FDC, D, all=TRUE)
D <- D[complete.cases(D),]

D <- SpatialPointsDataFrame(cbind(D$dec_long_va, D$dec_lat_va), D,
                            proj4string=LATLONG)
D <- spTransform(D, ALBEA)
XY <- coordinates(D)
D$east <- XY[,1]/1000; D$north <- XY[,2]/1000; rm(XY)

SO <- over(D, spDNI_1998to2009)
D$ANN_DNI <- SO$ANN_DNI
D$JAN <- SO$JAN; D$FEB <- SO$FEB
D$MAR <- SO$MAR; D$APR <- SO$APR
D$MAY <- SO$MAY; D$JUN <- SO$JUN
D$JUL <- SO$JUL; D$AUG <- SO$AUG
D$SEP <- SO$SEP; D$OCT <- SO$OCT
D$NOV <- SO$NOV; D$DEC <- SO$DEC
rm(SO)

bnd <- list(x=bnd_poly_aea[,1]/1000, y=bnd_poly_aea[,2]/1000)
bnd <- list(bnd)
plot(bnd[[1]]$x,bnd[[1]]$y,type="l", col=8, lwd=.6)
points(D$east, D$north, pch=4, lwd=.5, cex=0.9, col=8)

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
x <- knots$x; y <- knots$y
#knots <- data.frame(x=c(x, 1375), y=c(y, 600)); rm(x,y)
points(knots$x, knots$y, pch=16, cex=1.1, col=4)

duan_smearing_estimator <- function(model) { sum(10^residuals(model))/length(residuals(model)) }

z <- log10(D$L1)
x <- D$east; y <- D$north
L1   <- gam(z~te(CDA, ANN_DNI, bs="tp")+
              s(MAY, DEC, bs="tp", k=5)+
              s(x, y, bs="so", xt=list(bnd=bnd)), knots=knots,
              family="gaussian", data=D)
L1$duan_smearing <- duan_smearing_estimator(L1)
#pdf("L1.pdf", useDingbats=FALSE)
  plot(z, fitted.values(L1))
  abline(0,1)
  plot(L1, scheme=2)
  points(x, y, pch=4, lwd=.5, cex=0.9, col=8)
  points(knots$x, knots$y, pch=16, cex=1.1, col=4)
  text(100, 500, "L1")
#dev.off()
