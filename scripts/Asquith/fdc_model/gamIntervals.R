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
