
"ksvm2me" <- function(svm, obs, ...) {
  mysc <- slot(svm, "scaling")
  fit <- as.vector(mysc$y.scale$"scaled:scale" * predict(svm)+
                   mysc$y.scale$"scaled:center")
  n <- length(obs)
  if(length(fit) != n) {
     warning("predicted must be the same length as obs")
     return()
  }
  zz <- data.frame(obs=obs, fit=fit, svm.residuals=fit - obs)
  RSE <- sum((obs - fit)^2); RSE <- sqrt(RSE/length(obs))

  denom <- sum((obs - mean(obs))^2)
  if(denom != 0) {
     NSE <- 1 - sum((obs - fit)^2)/denom
  } else {
     NSE <- NA
     warning("'sum((obs - mean(obs))^2) == 0', not possible to compute 'NSE'")
  }
  ixin  <- slot(svm,"alphaindex")
  ixout <- (1:length(obs))[-ixin]
  zz$inSVM <- TRUE
  zz$inSVM[ixout] <- FALSE
  attr(zz, "RSE") <- RSE
  attr(zz, "NSE") <- NSE
  attr(zz, "frac.in" ) <- length(ixin) /n
  attr(zz, "frac.out") <- length(ixout)/n
  return(zz)
}



library(kernlab)
load("./Models.RData")
AD <- read.table("../sitefile/ExtraSiteFileStuff.txt", header=TRUE,
                 colClasses ="character")


Z <- D

Y <- log10(Z$L1)
sites <- D$site_no
area  <- D$basin_area
slope <- D$basin_slope
ppt   <- D$ppt_mean
temp  <- D$temp_mean
rad   <- D$dni_ann
dev   <- D$developed
flood  <- D$flood_storage
grass  <- D$grassland
decade <- as.numeric(D$decade)
n <- length(Y)
# RVMs are O(n^3) but result in radically fewer support vectors (prediction faster?)
RVM <- rvm(Y[1:n]~area[1:n]+slope[1:n]+ppt[1:n]+temp[1:n]+rad[1:n]+flood[1:n]+dev[1:n]+decade[1:n])
# SVMs are O(n^2) but results in more supports than RVMs, don't sweat RVM taking more time
SVM <- ksvm(Y[1:n]~area[1:n]+slope[1:n]+ppt[1:n]+temp[1:n]+rad[1:n]+flood[1:n]+dev[1:n]+decade[1:n],
            C=1, epsilon=.3)
SVMresults <- ksvm2me(SVM, Y)

plot(Y, SVMresults$fit, col=rgb(0,0,0.5,0.2))
points(Y, predict(RVM), col=rgb(0.5,0,0,0.2))

ix <- RVindex(RVM)
points(Y[ix], predict(RVM)[ix], col=2, pch=16)

points(Y[SVMresults$inSVM], SVMresults$fit[SVMresults$inSVM], pch=16, cex=0.5)

#https://cran.r-project.org/web/packages/kernlab/vignettes/kernlab.pdf
#Furthermore, unlike the SVM classifier, the non-zero weights in the RVM are not
# associated with examples close to the decision boundary, but rather appear to
# represent “prototypical” examples. These examples are termed relevance vectors


ix <- D$decade == "1950"; the.sites <- D$site_no[ix]
dRVM <- rvm(Y[ix]~area[ix]+slope[ix]+ppt[ix]+temp[ix]+rad[ix]+flood[ix])
ix <- RVindex(dRVM)
rsites <- the.sites[ix]
ix <- D$decade == "1960"; the.sites <- D$site_no[ix]
dRVM <- rvm(Y[ix]~area[ix]+slope[ix]+ppt[ix]+temp[ix]+rad[ix]+flood[ix])
ix <- RVindex(dRVM)
rsites <- c(rsites,the.sites[ix])
ix <- D$decade == "1970"; the.sites <- D$site_no[ix]
dRVM <- rvm(Y[ix]~area[ix]+slope[ix]+ppt[ix]+temp[ix]+rad[ix]+flood[ix])
ix <- RVindex(dRVM)
rsites <- c(rsites,the.sites[ix])
ix <- D$decade == "1980"; the.sites <- D$site_no[ix]
dRVM <- rvm(Y[ix]~area[ix]+slope[ix]+ppt[ix]+temp[ix]+rad[ix]+flood[ix])
ix <- RVindex(dRVM)
rsites <- c(rsites,the.sites[ix])
ix <- D$decade == "1990"; the.sites <- D$site_no[ix]
dRVM <- rvm(Y[ix]~area[ix]+slope[ix]+ppt[ix]+temp[ix]+rad[ix]+flood[ix])
ix <- RVindex(dRVM)
rsites <- c(rsites,the.sites[ix])
ix <- D$decade == "2000"; the.sites <- D$site_no[ix]
dRVM <- rvm(Y[ix]~area[ix]+slope[ix]+ppt[ix]+temp[ix]+rad[ix]+flood[ix])
ix <- RVindex(dRVM)
rsites <- c(rsites,the.sites[ix])
rsites <- unique(rsites)


ix <- D$decade == "1950"; the.sites <- D$site_no[ix]
dSVM <- ksvm(Y[ix]~area[ix]+slope[ix]+ppt[ix]+temp[ix]+rad[ix]+flood[ix])
ix <- ksvm2me(dSVM, Y[ix])$inSVM
ssites <- the.sites[ix]
ix <- D$decade == "1960"; the.sites <- D$site_no[ix]
dSVM <- ksvm(Y[ix]~area[ix]+slope[ix]+ppt[ix]+temp[ix]+rad[ix]+flood[ix])
ix <- ksvm2me(dSVM, Y[ix])$inSVM
ssites <- c(ssites,the.sites[ix])
ix <- D$decade == "1970"; the.sites <- D$site_no[ix]
dSVM <- ksvm(Y[ix]~area[ix]+slope[ix]+ppt[ix]+temp[ix]+rad[ix]+flood[ix])
ix <- ksvm2me(dSVM, Y[ix])$inSVM
ssites <- c(ssites,the.sites[ix])
ix <- D$decade == "1980"; the.sites <- D$site_no[ix]
dSVM <- ksvm(Y[ix]~area[ix]+slope[ix]+ppt[ix]+temp[ix]+rad[ix]+flood[ix])
ix <- ksvm2me(dSVM, Y[ix])$inSVM
ssites <- c(ssites,the.sites[ix])
ix <- D$decade == "1990"; the.sites <- D$site_no[ix]
dSVM <- ksvm(Y[ix]~area[ix]+slope[ix]+ppt[ix]+temp[ix]+rad[ix]+flood[ix])
ix <- ksvm2me(dSVM, Y[ix])$inSVM
ssites <- c(ssites,the.sites[ix])
ix <- D$decade == "2000"; the.sites <- D$site_no[ix]
dSVM <- ksvm(Y[ix]~area[ix]+slope[ix]+ppt[ix]+temp[ix]+rad[ix]+flood[ix])
ix <- ksvm2me(dSVM, Y[ix])$inSVM
ssites <- c(ssites,the.sites[ix])
ssites <- unique(ssites)


length(AD$site_no) # [1] 956
AD$active <- FALSE
AD$active[AD$end_cal_year >= 2016] <- TRUE
sum(AD$active) # [1] 673
# 673/956 thus 0.7039749 of the network is active


plot(D, lwd=0.4, cex=0)
for(site in AD$site_no) {
  active <- AD$active[AD$site_no == site]
  plot(D[D$site_no == site,], pch=2, cex=.5+0.5*as.numeric(active), col=2+as.numeric(active), lwd=0.7, add=TRUE)
}
for(site in rsites) {
#  plot(D[D$site_no == site,], pch=16, col=2, add=TRUE)
}
for(site in ssites) {
  active <- AD$active[AD$site_no == site]
  if(active) next
  plot(D[D$site_no == site,], pch=1, cex=1.5-as.numeric(active), col=4, lwd=0.7, add=TRUE)
}


Z <- D
Y <- as.factor(Z$nzero == 0)
dat <- data.frame(area=area,slope=slope,ppt=ppt,temp=temp,rad=rad,grass=grass,dev=dev,decade=decade)
SVM <- ksvm(Y~area+slope+ppt+temp+rad+grass+dev+decade, data=dat)
plot(Y,fitted(SVM))
sum(as.logical(Y) == as.logical(fitted(SVM)))/length(Y)
