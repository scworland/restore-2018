library(kernkab)

Y <- log10(D$L1)
sites <- D$site_no
area <- D$basin_area
slope <- D$basin_slope
ppt <- D$ppt_mean
temp <- D$temp_mean
rad <- D$dni_ann
dev <- D$developed
fs  <- D$flood_storage
decade <- as.numeric(D$decade)
n <- length(Y)
RVM <- rvm(Y[1:n]~A[1:n]+S[1:n]+P[1:n]+T[1:n]+R[1:n]+decade[1:n])

SVM <- ksvm(Y[1:n]~A[1:n]+S[1:n]+P[1:n]+T[1:n]+R[1:n]+decade[1:n], epsilon=0.2)


"ksvmHelper" <- function(svm, obs, ...) {
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
     warning("'sum((obs - mean(obs))^2)=0', not possible to compute 'NSE'")
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

SVMresults <- ksvmHelper(SVM, Y)

plot(Y, SVMresults$fit, col=rgb(0,0,0.5,0.2))
points(Y, predict(RVM), col=rgb(0.5,0,0,0.2))

ix <- RVindex(RVM)
points(Y[ix], predict(RVM)[ix], col=2, pch=16)

points(Y[SVMresults$inSVM], SVMresults$fit[SVMresults$inSVM])

#https://cran.r-project.org/web/packages/kernlab/vignettes/kernlab.pdf
#Furthermore, unlike the SVM classifier, the non-zero weights in the RVM are not associated
# with examples close to the decision boundary, but rather appear to represent “prototypical”
# examples. These examples are termed relevance vectors


ix <- D$decade == "1950"; the.sites <- D$site_no[ix]
RVM <- rvm(Y[ix]~area[ix]+slope[ix]+ppt[ix]+temp[ix]+rad[ix])
ix <- RVindex(RVM)
rsites <- the.sites[ix]
ix <- D$decade == "1960"; the.sites <- D$site_no[ix]
RVM <- rvm(Y[ix]~area[ix]+slope[ix]+ppt[ix]+temp[ix]+rad[ix])
ix <- RVindex(RVM)
rsites <- c(rsites,the.sites[ix])
ix <- D$decade == "1970"; the.sites <- D$site_no[ix]
RVM <- rvm(Y[ix]~area[ix]+slope[ix]+ppt[ix]+temp[ix]+rad[ix])
ix <- RVindex(RVM)
rsites <- c(rsites,the.sites[ix])
ix <- D$decade == "1980"; the.sites <- D$site_no[ix]
RVM <- rvm(Y[ix]~area[ix]+slope[ix]+ppt[ix]+temp[ix]+rad[ix])
ix <- RVindex(RVM)
rsites <- c(rsites,the.sites[ix])
ix <- D$decade == "1990"; the.sites <- D$site_no[ix]
RVM <- rvm(Y[ix]~area[ix]+slope[ix]+ppt[ix]+temp[ix]+rad[ix])
ix <- RVindex(RVM)
rsites <- c(rsites,the.sites[ix])
ix <- D$decade == "2000"; the.sites <- D$site_no[ix]
RVM <- rvm(Y[ix]~area[ix]+slope[ix]+ppt[ix]+temp[ix]+rad[ix])
ix <- RVindex(RVM)
rsites <- c(rsites,the.sites[ix])
rsites <- unique(rsites)

plot(D, lwd=0.5)
for(site in rsites) {
  plot(D[D$site_no == site,], pch=16, col=2, add=TRUE)
}

