
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
            C=1, epsilon=.2)
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
dSVM <- ksvm(Y[ix]~area[ix]+slope[ix]+ppt[ix]+temp[ix]+rad[ix]+flood[ix],C=1, epsilon=.2)
ix <- ksvm2me(dSVM, Y[ix])$inSVM
ssites <- the.sites[ix]
ix <- D$decade == "1960"; the.sites <- D$site_no[ix]
dSVM <- ksvm(Y[ix]~area[ix]+slope[ix]+ppt[ix]+temp[ix]+rad[ix]+flood[ix],C=1, epsilon=.2)
ix <- ksvm2me(dSVM, Y[ix])$inSVM
ssites <- c(ssites,the.sites[ix])
ix <- D$decade == "1970"; the.sites <- D$site_no[ix]
dSVM <- ksvm(Y[ix]~area[ix]+slope[ix]+ppt[ix]+temp[ix]+rad[ix]+flood[ix],C=1, epsilon=.2)
ix <- ksvm2me(dSVM, Y[ix])$inSVM
ssites <- c(ssites,the.sites[ix])
ix <- D$decade == "1980"; the.sites <- D$site_no[ix]
dSVM <- ksvm(Y[ix]~area[ix]+slope[ix]+ppt[ix]+temp[ix]+rad[ix]+flood[ix],C=1, epsilon=.2)
ix <- ksvm2me(dSVM, Y[ix])$inSVM
ssites <- c(ssites,the.sites[ix])
ix <- D$decade == "1990"; the.sites <- D$site_no[ix]
dSVM <- ksvm(Y[ix]~area[ix]+slope[ix]+ppt[ix]+temp[ix]+rad[ix]+flood[ix],C=1, epsilon=.2)
ix <- ksvm2me(dSVM, Y[ix])$inSVM
ssites <- c(ssites,the.sites[ix])
ix <- D$decade == "2000"; the.sites <- D$site_no[ix]
dSVM <- ksvm(Y[ix]~area[ix]+slope[ix]+ppt[ix]+temp[ix]+rad[ix]+flood[ix],C=1, epsilon=.2)
ix <- ksvm2me(dSVM, Y[ix])$inSVM
ssites <- c(ssites,the.sites[ix])
ssites <- unique(ssites)


length(AD$site_no) # [1] 956
AD$active <- FALSE
AD$active[AD$end_cal_year >= 2016] <- TRUE
sum(AD$active) # [1] 673
# 673/956 thus 0.7039749 of the network is active

pdf("PPLOfit_junk.pdf", useDingbats=FALSE, width=11, height=10)
  plot(spRESTORE_MGCV_BND)  # by creation of the PDF, we can get a handle on a global
  usr <- par()$usr # setting of the plotting limits by preserving the usr.
dev.off()
unlink("PPLOfit_junk.pdf")  # just quietly throw the file away

pdf("svm_rvm_inactivegages.pdf", useDingbats=FALSE, width=11, height=10)
map_base(xlim=usr[1:2], ylim=usr[3:4])
#plot(D, lwd=0.4, cex=0)
for(site in AD$site_no) {
  active <- AD$active[AD$site_no == site]
  plot(D[D$site_no == site,], pch=2+15*as.numeric(!active), cex=.7+0.2*as.numeric(! active), col=6-3*as.numeric(active), lwd=0.5, add=TRUE)
}
for(site in rsites) {
  active <- AD$active[AD$site_no == site]
  if(active) next
  plot(D[D$site_no == site,], pch=1, cex=1.5-as.numeric(active), col=2, lwd=0.7, add=TRUE)
}
for(site in ssites) {
  active <- AD$active[AD$site_no == site]
  if(active) next
  plot(D[D$site_no == site,], pch=1, cex=1.5-as.numeric(active), col=4, lwd=0.7, add=TRUE)
}
legend(-10000, 580000,
       c("Active (2016+) streamgage with at least one decade of record (1950s to 2000s)",
         "Inactive (2016-) streamgage with at least one decade of record (1950s to 2000s)",
         "Inactive site needed by RVM to estimate the log10(decadal mean nonzero streamflow)",
         "Inactive site needed by SVM to estimate the log10(decadal mean nonzero streamflow)"),
       cex=0.9, bty="n",
       pch=c(2,17,1,1), col=c(3,6,2,4))
map_annotation()
dev.off()
#system(paste0("pdfcrop --margins '-46 -110 -43 0' --clip svm_rvm_inactivegages.pdf svm_rvm_inactivegages.pdf"))

Z <- D
Y <- as.factor(Z$nzero == 0)
dat <- data.frame(area=area,slope=slope,ppt=ppt,temp=temp,rad=rad,grass=grass,dev=dev,decade=decade)
SVM <- ksvm(Y~area+slope+ppt+temp+rad+grass+dev+decade, data=dat,C=1, epsilon=.2)
plot(Y,fitted(SVM))
sum(as.logical(Y) == as.logical(fitted(SVM)))/length(Y)
