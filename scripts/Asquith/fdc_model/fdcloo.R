#source("fdcest.R")

library(mcparallelDo)


cvPPLOgo <- function(parent=PPLOdf, sigma=pplo.sigma, sites_to_fill=sites_to_fill) {
  sites <- decades <- values <- ests <- estlwrs <- estuprs <- NULL
  estft <- estlwrft <- estuprft <- NULL
  rse <- rsecv <- biascv <- 0; i <- 0
  for(site in unique(D$site_no)) { i <- i + 1
    Z <- D[D$site_no != site,]; x <- Z$x; y <- Z$y
    Z$left.threshold <-  log10(rep(0, length(Z$nzero)))
    Z$right.threshold <- log10(Z$n)
    Z$flowtime <- log10(Z$n - Z$nzero)
    model <- gam(flowtime~s(basin_area)+
                    s(ppt_mean, k=5)+s(temp_mean, k=4)+s(dni_ann, k=7)+
                    developed+grassland+
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
    message("cvPPLO: ",site, ", Bias=",biascv[i],", RSE=",rse[i],", RSEcv=",rsecv[i])
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
    model   <- gam(z~s(basin_area)+ s(basin_slope, k=5)+
               s(ppt_mean, k=5) + s(dni_ann, k=7)+
               developed+
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
    message("cvL1: ",site, ", Bias=",biascv[i],", RSE=",rse[i],", RSEcv=",rsecv[i])
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



cvT2go <- function(parent=T2df) {
  sites <- decades <- values <- ests <- estlwrs <- estuprs <- NULL
  rse <- rsecv <- biascv <- 0; i <- 0
  for(site in unique(D$site_no)) { i <- i + 1
    Z <- D[D$site_no != site,]; x <- Z$x; y <- Z$y
    z <- Z$L2/Z$L1 # LOOK HERE
    model   <- gam(z~s(basin_area) + s(basin_slope, k=5)+
               s(ppt_mean, k=5) + s(temp_mean, k=4) + s(dni_ann, k=7)+
               developed+
               s(flood_storage, k=5)+
               decade-1+
               s(x, y), knots=knots, data=Z, #, bs="so", xt=list(bnd=bnd)
               family="gaussian")
    rse[i] <- sqrt(summary(model)$scale)
    tmp <- D[D$site_no == site,]; val <- tmp$L2/tmp$L1 # LOOK HERE
    pgk <- predict(model, newdata=tmp, se.fit=TRUE)
    pgk <- gamIntervals(pgk, gam=model, interval="prediction")
    res <- (val - pgk$fit); biascv[i] <- mean(res)
    res <- sum(res^2); rsecv[i] <- sqrt(mean(res))
    message("cvT2: ",site, ", Bias=",biascv[i],", RSE=",rse[i],", RSEcv=",rsecv[i])
    if(is.null(values)) {
      sites <- tmp$site_no
      decades <- tmp$decade
      values <- val
      ests <- pgk$fit
      estlwrs <- pgk$lwr
      estuprs <- pgk$upr
    } else {
      s <- (length(values)+1); ss <- s:(s+length(val)-1)
      sites[ss] <- tmp$site_no
      decades[ss] <- tmp$decade
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
  zz <- data.frame(site_no_bak=sites, decade_bak=decades, orgfit=values,
             loo_est_lwr_T2=estlwrs, loo_est_T2=ests, loo_est_upr_T2=estuprs,
             stringsAsFactors = FALSE)
  return(zz)
}




cvT3go <- function(parent=T3df) {
  sites <- decades <- values <- ests <- estlwrs <- estuprs <- NULL
  rse <- rsecv <- biascv <- 0; i <- 0
  for(site in unique(D$site_no)) { i <- i + 1
    Z <- D[D$site_no != site,]; x <- Z$x; y <- Z$y
    z <- Z$T3 # LOOK HERE
    model   <- gam(z~s(basin_area) + s(basin_slope, k=5)+
               s(ppt_mean, k=5) + s(temp_mean, k=4) + s(dni_ann, k=7)+
               s(flood_storage, k=5)+
               decade-1+
               s(x, y), knots=knots, data=Z, #, bs="so", xt=list(bnd=bnd)
               family="gaussian")
    rse[i] <- sqrt(summary(model)$scale)
    tmp <- D[D$site_no == site,]; val <- tmp$T3 # LOOK HERE
    pgk <- predict(model, newdata=tmp, se.fit=TRUE)
    pgk <- gamIntervals(pgk, gam=model, interval="prediction")
    res <- (val - pgk$fit); biascv[i] <- mean(res)
    res <- sum(res^2); rsecv[i] <- sqrt(mean(res))
    message("cvT3: ",site, ", Bias=",biascv[i],", RSE=",rse[i],", RSEcv=",rsecv[i])
    if(is.null(values)) {
      sites <- tmp$site_no
      decades <- tmp$decade
      values <- val
      ests <- pgk$fit
      estlwrs <- pgk$lwr
      estuprs <- pgk$upr
    } else {
      s <- (length(values)+1); ss <- s:(s+length(val)-1)
      sites[ss] <- tmp$site_no
      decades[ss] <- tmp$decade
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
  zz <- data.frame(site_no_bak=sites, decade_bak=decades, orgfit=values,
             loo_est_lwr_T3=estlwrs, loo_est_T3=ests, loo_est_upr_T3=estuprs,
             stringsAsFactors = FALSE)
  return(zz)
}



cvT4go <- function(parent=T4df) {
  sites <- decades <- values <- ests <- estlwrs <- estuprs <- NULL
  rse <- rsecv <- biascv <- 0; i <- 0
  for(site in unique(D$site_no)) { i <- i + 1
    Z <- D[D$site_no != site,]; x <- Z$x; y <- Z$y
    z <- Z$T4 # LOOK HERE
    model   <- gam(z~s(basin_area) + s(basin_slope, k=5)+
               s(ppt_mean, k=5) + s(temp_mean, k=4) + s(dni_ann, k=7)+
               s(flood_storage, k=5)+
               decade-1+
               s(x, y), knots=knots, data=Z, #, bs="so", xt=list(bnd=bnd)
               family="gaussian")
    rse[i] <- sqrt(summary(model)$scale)
    tmp <- D[D$site_no == site,]; val <- tmp$T4 # LOOK HERE
    pgk <- predict(model, newdata=tmp, se.fit=TRUE)
    pgk <- gamIntervals(pgk, gam=model, interval="prediction")
    res <- (val - pgk$fit); biascv[i] <- mean(res)
    res <- sum(res^2); rsecv[i] <- sqrt(mean(res))
    message("cvT4: ",site, ", Bias=",biascv[i],", RSE=",rse[i],", RSEcv=",rsecv[i])
    if(is.null(values)) {
      sites <- tmp$site_no
      decades <- tmp$decade
      values <- val
      ests <- pgk$fit
      estlwrs <- pgk$lwr
      estuprs <- pgk$upr
    } else {
      s <- (length(values)+1); ss <- s:(s+length(val)-1)
      sites[ss] <- tmp$site_no
      decades[ss] <- tmp$decade
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
  zz <- data.frame(site_no_bak=sites, decade_bak=decades, orgfit=values,
             loo_est_lwr_T4=estlwrs, loo_est_T4=ests, loo_est_upr_T4=estuprs,
             stringsAsFactors = FALSE)
  return(zz)
}


cvT5go <- function(parent=T5df) {
  sites <- decades <- values <- ests <- estlwrs <- estuprs <- NULL
  rse <- rsecv <- biascv <- 0; i <- 0
  for(site in unique(D$site_no)) { i <- i + 1
    Z <- D[D$site_no != site,]; x <- Z$x; y <- Z$y
    z <- Z$T5 # LOOK HERE
    model   <- gam(z~s(basin_area) + s(basin_slope, k=5)+
               s(temp_mean, k=4) + s(dni_ann, k=7)+
               s(flood_storage, k=5)+
               decade-1+
               s(x, y), knots=knots, data=Z, #, bs="so", xt=list(bnd=bnd)
               family="gaussian")
    rse[i] <- sqrt(summary(model)$scale)
    tmp <- D[D$site_no == site,]; val <- tmp$T5 # LOOK HERE
    pgk <- predict(model, newdata=tmp, se.fit=TRUE)
    pgk <- gamIntervals(pgk, gam=model, interval="prediction")
    res <- (val - pgk$fit); biascv[i] <- mean(res)
    res <- sum(res^2); rsecv[i] <- sqrt(mean(res))
    message("cvT5: ",site, ", Bias=",biascv[i],", RSE=",rse[i],", RSEcv=",rsecv[i])
    if(is.null(values)) {
      sites <- tmp$site_no
      decades <- tmp$decade
      values <- val
      ests <- pgk$fit
      estlwrs <- pgk$lwr
      estuprs <- pgk$upr
    } else {
      s <- (length(values)+1); ss <- s:(s+length(val)-1)
      sites[ss] <- tmp$site_no
      decades[ss] <- tmp$decade
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
  zz <- data.frame(site_no_bak=sites, decade_bak=decades, orgfit=values,
             loo_est_lwr_T5=estlwrs, loo_est_T5=ests, loo_est_upr_T5=estuprs,
             stringsAsFactors = FALSE)
  return(zz)
}




mcparallelDo({
  cvPPLO <- cvPPLOo <- cvPPLOgo()

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
  PPLOdf_loo$key <- PPLOdf$key <- NULL
  PPLOdf_loo$site_no_bak <- NULL
  PPLOdf_loo$decade_bak  <- NULL
  PPLOdf_loo$orgfit      <- NULL
  PPLOdf_loo
}, "looPPLOmodel")


mcparallelDo({
cvL1 <- cvL1o <- cvL1go()

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
L1df_loo$key <- L1df$key <- NULL
L1df_loo$site_no_bak <- NULL
L1df_loo$decade_bak  <- NULL
L1df_loo$orgfit      <- NULL
L1df_loo
}, "looL1model")


mcparallelDo({
  cvT2 <- cvT2o <- cvT2go()

  for(site in sites_to_fill) {
     tmp <- DDo[DDo$site_no == site,]
     df <- data.frame(site_no_bak=tmp$site_no, decade_bak=tmp$decade,
                      orgfit=tmp$L2/tmp$L1,
                      loo_est_lwr_T2=NA, loo_est_T2=NA, loo_est_upr_T2=NA,
                      stringsAsFactors=FALSE)
     cvT2 <- rbind(cvT2, df)
  }
  cvT2 <- cvT2[order(cvT2$site_no_bak, cvT2$decade_bak),]
  cvT2$key <-  paste(cvT2$site_no_bak, cvT2$decade_bak, sep=":")

  T2df$key <- paste(T2df$site_no, T2df$decade, sep=":")
  T2df_loo <- merge(T2df, cvT2, add=TRUE)
  T2df_loo$key <- T2df$key <- NULL
  T2df_loo$site_no_bak <- NULL
  T2df_loo$decade_bak  <- NULL
  T2df_loo$orgfit      <- NULL
  T2df_loo
}, "looT2model")

mcparallelDo({
 cvT3 <- cvT3o <- cvT3go()

  for(site in sites_to_fill) {
     tmp <- DDo[DDo$site_no == site,]
     df <- data.frame(site_no_bak=tmp$site_no, decade_bak=tmp$decade,
                      orgfit=tmp$T3,
                      loo_est_lwr_T3=NA, loo_est_T3=NA, loo_est_upr_T3=NA,
                      stringsAsFactors=FALSE)
     cvT3 <- rbind(cvT3, df)
  }
  cvT3 <- cvT3[order(cvT3$site_no_bak, cvT3$decade_bak),]
  cvT3$key <-  paste(cvT3$site_no_bak, cvT3$decade_bak, sep=":")

  T3df$key <- paste(T3df$site_no, T3df$decade, sep=":")
  T3df_loo <- merge(T3df, cvT3, add=TRUE)
  T3df_loo$key <- T3df$key <- NULL
  T3df_loo$site_no_bak <- NULL
  T3df_loo$decade_bak  <- NULL
  T3df_loo$orgfit      <- NULL
  T3df_loo
}, "looT3model")



mcparallelDo({
  cvT4 <- cvT4o <- cvT4go()

  for(site in sites_to_fill) {
     tmp <- DDo[DDo$site_no == site,]
     df <- data.frame(site_no_bak=tmp$site_no, decade_bak=tmp$decade,
                      orgfit=tmp$T4,
                      loo_est_lwr_T4=NA, loo_est_T4=NA, loo_est_upr_T4=NA,
                      stringsAsFactors=FALSE)
     cvT4 <- rbind(cvT4, df)
  }
  cvT4 <- cvT4[order(cvT4$site_no_bak, cvT4$decade_bak),]
  cvT4$key <-  paste(cvT4$site_no_bak, cvT4$decade_bak, sep=":")

  T4df$key <- paste(T4df$site_no, T4df$decade, sep=":")
  T4df_loo <- merge(T4df, cvT4, add=TRUE)
  T4df_loo$key <- T4df$key <- NULL
  T4df_loo$site_no_bak <- NULL
  T4df_loo$decade_bak  <- NULL
  T4df_loo$orgfit      <- NULL
 T4df_loo
}, "looT4model")


mcparallelDo({
  cvT5 <- cvT5o <- cvT5go()

  for(site in sites_to_fill) {
     tmp <- DDo[DDo$site_no == site,]
     df <- data.frame(site_no_bak=tmp$site_no, decade_bak=tmp$decade,
                      orgfit=tmp$T5,
                      loo_est_lwr_T5=NA, loo_est_T5=NA, loo_est_upr_T5=NA,
                      stringsAsFactors=FALSE)
     cvT5 <- rbind(cvT5, df)
  }
  cvT5 <- cvT5[order(cvT5$site_no_bak, cvT5$decade_bak),]
  cvT5$key <-  paste(cvT5$site_no_bak, cvT5$decade_bak, sep=":")

  T5df$key <- paste(T5df$site_no, T5df$decade, sep=":")
  T5df_loo <- merge(T5df, cvT5, add=TRUE)
  T5df_loo$key <- T5df$key <- NULL
  T5df_loo$site_no_bak <- NULL
  T5df_loo$decade_bak  <- NULL
  T5df_loo$orgfit      <- NULL
  T5df_loo
}, "looT5model")


while(1) {
   Sys.sleep(10)
   me <- mcparallelDoCheck()
   print(me)
   if(length(me) == 0) break
   Sys.sleep(50)
}

looOverL1model <-
  data.frame(comid=looL1model$comid, site_no=looL1model$site_no,
             huc12=looL1model$huc12, decade=looL1model$decade,
             dec_long_va=looL1model$dec_long_va,
             dec_lat_va=looL1model$dec_lat_va,
             in_model_pplo=looPPLOmodel$in_model_pplo,
             in_model_L1=looL1model$in_model_L1,
             overL1=(1-looPPLOmodel$pplo)*looL1model$L1,
             stringsAsFactors=FALSE)
looOverL1model$est_lwr_overL1 <- (1-looPPLOmodel$est_upr_pplo)*looL1model$est_lwr_L1*looL1model$bias_corr
looOverL1model$est_overL1     <- (1-looPPLOmodel$est_pplo)    *looL1model$est_L1*looL1model$bias_corr
looOverL1model$est_upr_overL1 <- (1-looPPLOmodel$est_lwr_pplo)*looL1model$est_upr_L1*looL1model$bias_corr
looOverL1model$loo_est_lwr_overL1 <- (1-looPPLOmodel$loo_est_upr_pplo)*looL1model$loo_est_lwr_L1*looL1model$bias_corr
looOverL1model$loo_est_overL1     <- (1-looPPLOmodel$loo_est_pplo)    *looL1model$loo_est_L1*looL1model$bias_corr
looOverL1model$loo_est_upr_overL1 <- (1-looPPLOmodel$loo_est_lwr_pplo)*looL1model$loo_est_upr_L1*looL1model$bias_corr

assign("last.warning", NULL, envir = baseenv())
write_feather(looOverL1model, "all_gage_looest_overL1.feather") # must then run fdcloo_post.R
write_feather(looPPLOmodel, "all_gage_looest_pplo.feather")
write_feather(looL1model,   "all_gage_looest_L1.feather"  )
write_feather(looT2model,   "all_gage_looest_T2.feather"  )
write_feather(looT3model,   "all_gage_looest_T3.feather"  )
write_feather(looT4model,   "all_gage_looest_T4.feather"  )
write_feather(looT5model,   "all_gage_looest_T5.feather"  )
