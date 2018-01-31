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
knots <- knots[-c(1, 3, 5, 6, 13, 10, 11, 12, 16, 21, 22, 27, 28, 29),]
x <- knots$x; y <- knots$y
knots <- data.frame(x=c(x, 1375), y=c(y, 600)); rm(x,y)
points(knots$x, knots$y, pch=16, cex=1.1, col=4)

knots_pplo <- knots[-17, ]
knots_pplo <- knots
# The as.numeric() is needed head of pplo use later
plogit <- function(eta) as.numeric(exp(eta)/(exp(eta)+1))
duan_smearing_estimator <- function(model) { sum(10^residuals(model))/length(residuals(model)) }

DD$x <- DD$east; DD$y <- DD$north

#dd <- slot(DD, name="data")
#for(i in c(4:5, 45:79))
#n     <- aggregate(DD$n,     by=list(DD$site_no), sum)$x
#nzero <- aggregate(DD$nzero, by=list(DD$site_no), sum)$x










DD$flood_storage <- (DD$tot_nid_storage - DD$tot_norm_storage)/10^DD$CDA
DD[DD$flood_storage < 0,] # two sites: 02295420 and 02296750
DD$flood_storage[DD$flood_storage < 0] <- (DD$tot_norm_storage[DD$flood_storage < 0] -
                                            DD$tot_nid_storage[DD$flood_storage < 0]) /
                                                     10^DD$CDA[DD$flood_storage < 0]
DD$flood_storage[DD$flood_storage > 3000] <- 3000
DD$flood_storage <- log10(DD$flood_storage+0.01)

DD$tot_basin_area  <- log10(DD$tot_basin_area)
DD$tot_basin_slope <- log10(DD$tot_basin_slope)
sites_of_area_bust <- unique(DD$site_no[abs(DD$tot_basin_area - DD$CDA) > 1/3])
for(site in sites_of_area_bust) {
  DD <- DD[DD$site_no != site,]
}

DD$decade <- as.factor(DD$decade)
DD$bedperm <- as.factor(DD$bedperm)
DD$aquifers <- as.factor(DD$aquifers)
DD$soller <- as.factor(DD$soller)
DD$hlr <- as.factor(DD$hlr)
DD$ecol3 <- as.factor(DD$ecol3)
DD$physio <- as.factor(DD$physio)
DD$statsgo <- as.factor(DD$statsgo)

DD$isWest <- 0
DD$isWest[DD$x < 0] <- 1
DD$isFL <- 0
DD$isFL[DD$x > 1200 & DD$y < 750] <- 1
DD$isWest <- as.logical(DD$isWest)
DD$isFL <- as.logical(DD$isFL)




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
DD$alt_physio[DD$physio == "cat_physio_12"] <- "Central Lowland"
DD$alt_physio[DD$physio == "cat_physio_13"] <- "Great Plains"
DD$alt_physio <- as.factor(DD$alt_physio)
DD$alt_physio <- relevel(DD$alt_physio, "other")





DD$texas_special <- 0
DD$texas_special[DD$site_no == "08156800"] <- 1
DD$texas_special[DD$site_no == "08155300"] <- 1
DD$texas_special[DD$site_no == "08155400"] <- 1
DD$texas_special[DD$site_no == "08181400"] <- 1
DD$texas_special[DD$site_no == "08184000"] <- 1
DD$texas_special[DD$site_no == "08185000"] <- 1
DD$texas_special[DD$site_no == "08190500"] <- 1
#DD$texas_special[DD$site_no == "08194200"] <- 1
DD$texas_special[DD$site_no == "08192000"] <- 1
DD$texas_special[DD$site_no == "08197500"] <- 1
DD$texas_special[DD$site_no == "08198500"] <- 1
DD$texas_special[DD$site_no == "08200700"] <- 1
DD$texas_special[DD$site_no == "08202700"] <- 1
#DD$texas_special[DD$site_no == "08205500"] <- 1
#DD$texas_special[DD$site_no == "08212400"] <- 1
#DD$texas_special[DD$site_no == "08205500"] <- 1
DD$texas_special <- as.factor(DD$texas_special)

D <- DD;

D <- D[D$site_no != "08066191", ] # Livingston Res Outflow Weir nr Goodrich, TX
D <- D[D$site_no != "08154510", ] # Colorado Rv bl Mansfield Dam, Austin, TX
D <- D[D$site_no != "08158000", ] # Colorado Rv at Austin, TX
D <- D[D$site_no != "08169000", ] # Comal River (major spring)
D <- D[D$site_no != "08170500", ] # San Marcos River (major spring)
D <- D[D$texas_special != "1",]


save(D, DD, knots, knots_pplo, bnd, file="DEMO.RData")

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
plot(bnd[[1]]$x,bnd[[1]]$y,type="l", col=8, lwd=.6)
points(D$east[D$nzero == 0], D$north[D$nzero == 0], pch=4, lwd=.5, cex=0.9, col=rgb(0,.5,0.5,.2))
points(D$east[D$nzero > 0 & D$nzero <= 800], D$north[D$nzero > 0 & D$nzero <= 800], pch=4, lwd=.5, cex=0.9, col=rgb(1,0,0.5,.2))
points(D$east[D$nzero > 2000], D$north[D$nzero > 2000], pch=16, lwd=.5, cex=0.9, col=rgb(0.5,0,1,.4))

Z <- D
Z <- Z[Z$nzero > 0,]
#x <- Z$east; y <- Z$north # Please see the extraction and setting to x,y
#CDA <- Z$CDA; ANN_DNI <- Z$ANN_DNI; MAY <- Z$MAY; DEC <- Z$DEC
Z$z <- log10(Z$nzero)
kap <- parkap(lmoms(Z$z))
plot(pp(Z$z), sort(Z$z)); FF <- nonexceeds(); lines(FF, qlmomco(FF, kap), col=2)
Z$z <- qnorm(plmomco(Z$z, kap))
#Z$z <- cbind(Z$n-Z$nzero, Z$n)

# I(2*asin(sqrt(developed/100)))+I(sqrt(tot_hdens))
#s(ANN_DNI, bs="cr", k=5)
PPLO <- gam(z~s(CDA, bs="cr", k=4)+
              s(ppt_mean)+alt_physio,
              knots=knots_pplo,
              data=Z)
PPLO <- gam(z~CDA+log10(ppt_mean)+alt_physio+decade,
              knots=knots_pplo,
              data=Z)


family <- "weibull"
CCC <- NULL
decades <- levels(D$decade)
for(i in 1:length(decades)) {
   decade <- decades[i]
   Z <-  D[D$decade == decade,]
   Z <- D
   cda <- Z$CDA
   developed <- 2*asin(sqrt(Z$developed/100))
   tot_hdens <-  log10(Z$tot_hdens)
   cultivated_cropland <- 2*asin(sqrt(Z$cultivated_cropland/100))
   ppt_mean <- log10(Z$ppt_mean)
   flood_storage <- Z$flood_storage; tot_basin_slope <- Z$tot_basin_slope
   Zc <- Surv(Z$n - Z$nzero, Z$nzero != 0, type="right")
   plot(survfit(Zc~1)); mtext(paste0("Decade ",decade))

   SM <- survreg(Zc~cda+ppt_mean+Z$ANN_DNI+tot_basin_slope+developed+cultivated_cropland+Z$decade, dist=family)
   sSM <- summary(SM)
   #print(sSM)

   PPo <- predict(SM)
   PP <- 10^PPo
   daysUntilzero <- Z$n - Z$nzero; nmax <- max(daysUntilzero)
   PP[PP > nmax] <- nmax
   PPc <- PP[PP != nmax]
   residuals <- daysUntilzero - PP
   lim <- range(c(daysUntilzero, PP))
   plot(daysUntilzero[PP != daysUntilzero], PP[PP != daysUntilzero],
        pch=16, col=rgb(0,0,0,.2), log="xy", xlim=lim, ylim=lim)
   abline(0,1); mtext(paste0("Decade ",decade))
   res <- residuals(SM)
   xlim <- range(PPo)
   plot(PPo, res, xlim=xlim, ylim=c(-1.5,0.5), type="n", xaxs="i", yaxs="i",
        xlab="log10(predicted days until zero flow)", ylab="Residual, log10")
   lines(rep(log10(3653), 2), par()$usr[3:4], lty=1, lwd=0.5)
   X <- 10^c(seq(par()$usr[1],par()$usr[2], by=.02), par()$usr[2])
   lines(log10(X), log10(X) - log10(X- 7*10), lty=4, lwd=0.5)
   lines(log10(X), log10(X) - log10(X+ 7*10), lty=4, lwd=0.5)
   #lines(log10(X), log10(X) - log10(X-14*10), lty=1, lwd=0.5)
   #lines(log10(X), log10(X) - log10(X+14*10), lty=1, lwd=0.5)
   lines(log10(X), log10(X) - log10(X-30*10), lty=1, lwd=0.5)
   lines(log10(X), log10(X) - log10(X+30*10), lty=1, lwd=0.5)
   lines(log10(X), log10(X) - log10(X-60*10), lty=2, lwd=0.5)
   lines(log10(X), log10(X) - log10(X+60*10), lty=2, lwd=0.5)
   tmp <- sapply(1:length(PPo), function(i) {
     alive <- ! as.logical(as.numeric(Zc[i])[2])
     col <- ifelse(PP[i] == daysUntilzero[i], rgb(0,0,1,.3), rgb(1,0,0,.3))
     if(alive) segments(x0=PPo[i], x1=PPo[i], y0=res[i], y1=0.5, col=col, lwd=0.8)
   })
   points(PPo[PP == daysUntilzero], res[PP == daysUntilzero], lwd=0.6, cex=0.6, col=4)
   points(PPo[PP != daysUntilzero], res[PP != daysUntilzero], lwd=0.6, cex=0.8, col=2)
   mtext(paste0("Decade ",decade))
   correct_decision_rate <- sum(PP == daysUntilzero)/length(PPo)

   CC <- as.data.frame(t(coefficients(SM)))
   tmp <- names(CC); tmp[1] <- "Intercept"; names(CC) <- tmp
   CC$n <- length(PPo); CC$correct_decision_rate <- correct_decision_rate
   CC$MAD <- mean(abs(residuals))
   CC$MAD_nonzero <- 10^mean(log10(abs(residuals[residuals != 0])))
   CC$AIC <- AIC(SM)
   print(CC)
   if(is.null(CCC)) {
      CCC <- CC
   } else {
      CCC <- rbind(CCC,CC)
   }
}
CCC$decade <- as.numeric(decades)
for(i in 1:length(CCC[1,])) {
  name <- names(CCC)[i]
  vals <- CCC[,i]
  wgts <- CCC$n
  message(name, "  coef=",weighted.mean(vals, wgts))
}







smTillZero <- function(cda, ppt_mean, may, dec, developed, region) {
   developed <- 2*asin(sqrt(developed/100))
   return(1.32594940801398 +0.0684549053731752*log10(cda)+
          0.275911326088305*log10(ppt_mean)+
          0.336112889429094*may+0.219503899716695*dec-
          0.0637326522016246*may*dec+0.0605904049010024*developed)
   reg <- sapply(region, function(r) {
          ifelse(r == "-27", -0.11031,
          ifelse(r == "-31", -0.23108,
          ifelse(r == "-34", -0.07966,
          ifelse(r == "-75", -0.16082, 0)))) })
   b <- -1.60021 + reg
   cda <- log10(cda)
   ppt_mean <- log10(ppt_mean)
   tmp <- 0.12726*cda + 0.71234*ppt_mean +0.64845*may
   tmp <- tmp + 0.50603*dec + 0.00226*developed
   tmp <- tmp + -0.12340*may*dec + b
   names(tmp) <- NULL
   return(tmp)
}

GG <- 10^predict(SM)
#GG[GG > 3653] <- 3653

G <- smTillZero(10^Z$CDA, Z$ppt_mean, Z$MAY, Z$DEC, Z$developed, Z$alt_ecol3)
G <- 10^G
#G[G > 3653] <- 3653
plot(G,GG)

G <- smTillZero(10^D$CDA, D$ppt_mean, D$MAY, D$DEC, D$developed, D$alt_ecol3)
G <- 10^G

plot(D$n - D$nzero, G, log="xy")
abline(0,1)







L1s <- T2s <- T3s <- T4s <- Ns <- coeBs <- coeCDAs <- ppts <- units <- coedevs <- rep(NA, length(levels(Z$decade)))
PPLOSenv <- new.env(); decades <- levels(Z$decade)
for(i in 1:length(decades)) {
  decade <- decades[i]
  file <- paste0("PPLO",decade,".pdf"); message(file)
  pdf(file, useDingbats=FALSE)
  Z <-  D[D$decade == decade,]; nZ <- Z[Z$nzero > 0,]
  n <- length(nZ$pplo)
  nZ$z <- log10(nZ$nzero); lmr <- lmoms(nZ$z)
  Ns[i] <- n; L1s[i] <- lmr$lambdas[1]
  T2s[i] <- lmr$ratios[2]; T3s[i] <- lmr$ratios[3]; T4s[i] <- lmr$ratios[4]
  kap <- parkap(lmr); print(kap$para)
  plot(pp(nZ$z), sort(nZ$z)); FF <- nonexceeds(); lines(FF, qlmomco(FF, kap), col=2)
  nZ$z <- plmomco(nZ$z, kap)
  nZ$z[nZ$z == 1] <- 0.9999; nZ$z[nZ$z == 0] <- 1-0.9999
  nZ$z <- qnorm(nZ$z)
  mtext(paste0("Decade ",decade))


  z <- nZ$z; CDA <- nZ$CDA; ppt_mean <- log10(nZ$ppt_mean); unitL1 <- log10(nZ$L1)
  ANN_DNI <- nZ$ANN_DNI; x <- nZ$x; y <- nZ$y; T3 <- nZ$T3
  developed <- 2*asin(sqrt(nZ$developed/100))
  PPLO <- gam(z~CDA+T3+s(ppt_mean)+developed+s(ANN_DNI)+s(x, y, bs="so", xt=list(bnd=bnd)), knots=knots)#++I(2*asin(sqrt(developed/100)))+alt_physio, data=nZ)
  PPLO <- gam(z~CDA+developed+
                s(ANN_DNI, bs="cr", k=8))
  coes <- coefficients(PPLO)[1:5]
  coeBs[i] <- coes[1]; coeCDAs[i] <- coes[2]; ppts[i] <- coes[3]; units[i] <- coes[4]; coedevs[i] <- coes[5]
  PPLO$the.Z <- nZ; PPLO$the.kap <- kap; PPLO$the.lmr <- lmr
  assign(decade, PPLO, envir=PPLOSenv)
  print(summary(PPLO))
  #plot(nZ$nzero, 10^qlmomco(pnorm(fitted.values(PPLO)), kap), log="xy"); abline(0,1)
  #mtext(paste0("Decade ",decade))
  plot(nZ$CDA, nZ$nzero, pch=16, col=8, cex=developed+0.4)
  points(nZ$CDA, 10^qlmomco(pnorm(fitted.values(PPLO)), kap), cex=developed+0.4)
  mtext(paste0("Decade ",decade))
  #plot(PPLO, scheme=2, residuals=TRUE, pch=16, cex=0.5)
  #points(nZ$x, nZ$y, pch=4, lwd=.5, cex=0.9, col=8)
  #points(knots$x, knots$y, pch=16, cex=1.1, col=4)
  #text(100, 500, "Z-score of number of noflow days", pos=4)
  dev.off()
}
DF <- data.frame(decade=decades, n=Ns, L1=L1s, T2=T2s, T3=T3s, T4=T4s,
                 intercept=coeBs, cda=coeCDAs, ppt_mean=ppts, unitL1CDA=units, dev=coedevs)
assign("data.summary", DF, envir=PPLOSenv)
save(PPLOSenv, file="PPLOS.RData")

weighted.mean(DF$L1, DF$n)
weighted.mean(DF$T2, DF$n)
weighted.mean(DF$T3, DF$n)
weighted.mean(DF$T4, DF$n)
weighted.mean(DF$intercept, DF$n)
weighted.mean(DF$cda, DF$n)
weighted.mean(DF$ppt_mean, DF$n)
weighted.mean(DF$unitL1CDA, DF$n)
weighted.mean(DF$dev, DF$n)

solve.pplo <- function(cda=NA, ppt_mean=NA, L1=NA, developed=NA) {
  intercept   <- -7.7014393
  cda.c       <-  0.6338891
  ppt_mean.c  <-  1.9915757
  unitL1CDA.c <- -1.0132245
  developed.c <-  0.1665043
  developed <- 2*asin(sqrt(developed/100))
  z.score <- intercept +  cda.c      * log10(cda)    + ppt_mean.c  * log10(ppt_mean) +
                          unitL1CDA.c* log10(L1/cda) + developed.c * developed
  lmr  <- lmomco::vec2lmom(c(2.256159, 0.2025797, -0.1334176, 0.08901554), lscale=FALSE)
  kap  <- lmomco::parkap(lmr)
  pplo <- lmomco::qlmomco(pnorm(z.score), kap)
  if(pplo < 0) pplo <- 0
  if(pplo > 365.25*10) {
     warning("it nevers flows, later computations will fall apart")
     pplo <- 1
  }
  return(10^pplo)
}
solve.pplo(cda=10^2.58, ppt_mean=1293.2, L1=263.45, developed=8.84)





  #plot(asin(sqrt(Z$pplo)), asin(sqrt(fitted.values(PPLO)))); abline(0,1)
  plot(Z$nzero, 10^qlmomco(pnorm(fitted.values(PPLO)), kap), log="xy"); abline(0,1)
  plot(Z$CDA, Z$pplo, pch=16, col=8)
  points(Z$CDA, fitted.values(PPLO))
  plot(PPLO, scheme=2)
  points(Z$x, Z$y, pch=4, lwd=.5, cex=0.9, col=8)
  points(knots_pplo$x, knots_pplo$y, pch=16, cex=1.1, col=4)
  text(100, 500, "Probability level of noflows", pos=4)
dev.off()

ZEROSenv <- new.env()
for(decade in levels(D$decade)) {
  file <- paste0("ZERO",decade,".pdf"); message(file)
  pdf(file, useDingbats=FALSE)
  Z <- D[D$decade == decade,]
  Z$z <- Z$prob_has_noflow
  ZERO <- gam(z~CDA+I(2*asin(sqrt(developed/100)))+
                s(ANN_DNI, bs="cr", k=6)+s(MAY, DEC, bs="tp", k=6),
              knots=knots_pplo, data=Z, family="quasibinomial")
  ZERO$the.Z <- Z
  assign(decade, ZERO, envir=ZEROSenv)
  print(ZERO$coefficients[1:3])

  fv <- fitted.values(ZERO); atem <- rep(0.04, length(Z$z)); atem[fv > 0.5] <- 0.96
  plot(Z$prob_has_noflow, fv, xlab="Has some no flow", ylab="Fitted value")
  mtext(paste0("Decade ",decade))
  plot(Z$ANN_DNI, jitter(Z$z, amount=0.02), pch=16, col=rgb(0,0,1,.3),
       xlab="Annual irradiance (not actually in the model)",
       ylab="Has some no flow (red is the model inset a few percent)")
  n <- length(Z$prob_has_noflow)
  j <- sum(Z$prob_has_noflow[Z$prob_has_noflow == 0 & atem == 0.04]+1)
  k <- sum(Z$prob_has_noflow[Z$prob_has_noflow == 1 & atem == 0.96])
  coherence <- 100*(j + k)/n; message("coherence ", as.integer(coherence))
  mtext(paste0("Decade ",decade, " with coeherence = ", round(coherence, digits=2),
               " percent (correct decision)"))
  points(Z$ANN_DNI, jitter(atem, amount=0.02), pch=16, col=rgb(1,0,0,.3))

  plot(ZERO, scheme=2, residuals=TRUE, pch=16, cex=0.5)
  #points(Z$x, Z$y, pch=4, lwd=.5, cex=0.9, col=8)
  #points(knots_pplo$x, knots_pplo$y, pch=16, cex=1.1, col=4)
  text(100, 550, decade)
  text(100, 500, "Probability model of site having at least one noflow", pos=4, cex=0.7)
  dev.off()
}



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
