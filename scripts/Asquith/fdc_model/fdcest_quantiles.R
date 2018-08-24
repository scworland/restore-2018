



z <- log10(D$f50+1)      # ---------------------------
F50   <- gam(z~s(basin_area) +
               s(ppt_mean, k=5) + s(temp_mean, k=4) + s(dni_ann, k=7)+
               developed+
               mixed_forest+shrubland+
               decade-1+
               s(x, y), knots=knots, data=D, #, bs="so", xt=list(bnd=bnd)
             family="gaussian")
pdf("F50.pdf", useDingbats=FALSE)
plot(z, fitted.values(F50))
abline(0,1)
plot(F50, scheme=2, residuals=TRUE)
points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
points(knots$x, knots$y, pch=16, cex=1.1, col=4)
text(100, 500, "F50")
dev.off()



PGAM <- gamIntervals(predict(F50, se.fit=TRUE), gam=F50, interval="prediction")
#plot(F50$hat, (PGAM$se.fit/PGAM$residual.scale)^2)
#abline(0,1)
sigma <- PGAM$residual.scale[1]
# Terms invert in upper/lower meaning and hence the flipping during data.frame construction.
F50df <- data.frame(comid=Z$comid, site_no=Z$site_no, huc12=Z$huc12,
                     decade=Z$decade,
                     dec_long_va=Z$dec_long_va, dec_lat_va=Z$dec_lat_va,
                     in_model_f50=TRUE, f50=Z$f50,
                     est_lwr_f50=10^PGAM$lwr-1,
                     est_f50    =10^PGAM$fit-1,
                     est_upr_f50=10^PGAM$upr-1, stringsAsFactors=FALSE)
F50df$rse_f50 <- sigma
F50df$se.fit_f50 <- PGAM$se.fit

sites_to_fill <- unique(c(sites_of_area_bust, DDo$site_no[DDo$ed_rch_zone == 1]))
i <- 0
for(site in sites_to_fill) {
  for(decade in as.character(DDo$decade[DDo$site_no == site])) {
    i <- i + 1
    message(site," ", decade, " ", i)
    tmp <- data.frame(comid=DDo$comid[DDo$site_no == site & DDo$decade == decade],
                      site_no=site,
                      huc12=DDo$huc12[DDo$site_no == site & DDo$decade == decade],
                      decade=decade,
                      dec_long_va=DDo$dec_long_va[DDo$site_no == site & DDo$decade == decade],
                      dec_lat_va=DDo$dec_lat_va[DDo$site_no == site & DDo$decade == decade],
                      in_model_f50=FALSE,
                      f50=DDo$f50[DDo$site_no == site & DDo$decade == decade],
                      est_lwr_f50=NA, est_f50=NA, est_upr_f50=NA,
                      rse_f50=NA, se.fit_f50=NA, stringsAsFactors=FALSE)
    F50df <- rbind(F50df, tmp)
  }
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(F50, newdata=tmp, se.fit=TRUE)
}
F50df <- F50df[order(F50df$site_no, F50df$decade),]

for(site in sites_to_fill) {
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(F50, newdata=tmp, se.fit=TRUE)
  pgk <- gamIntervals(jnk, gam=F50, interval="prediction")
  df <- data.frame(comid=tmp$comid, site_no=tmp$site_no, huc12=tmp$huc12,
                   decade=tmp$decade,
                   dec_long_va=tmp$dec_long_va, dec_lat_va=tmp$dec_lat_va,
                   in_model_f50=0,
                   f50=tmp$f50,
                   est_lwr_f50=10^pgk$lwr-1,
                   est_f50    =10^pgk$fit-1,
                   est_upr_f50=10^pgk$upr-1, stringsAsFactors=FALSE)
  df$rse_f50 <- sigma
  df$se.fit_f50 <- pgk$se.fit
  F50df[F50df$site_no == site,] <- df
}
F50df$est_f50[F50df$est_f50 < 0] <- 0
F50df$est_lwr_f50[F50df$est_lwr_f50 < 0] <- 0
F50df$est_upr_f50[F50df$est_upr_f50 < 0] <- 0

#write_feather(F50df, "all_gage_est_f50.feather")



#sink("right_tail_flowing_fdc.txt")
z <- log10(D$f90+1)      # ---------------------------
F90   <- gam(z~s(basin_area)  + s(flood_storage, k=7) +
              s(ppt_mean, k=5) + s(dni_ann, k=7)+
              developed+
              mixed_forest+shrubland+
              decade-1+
              s(x, y), knots=knots, data=D, #, bs="so", xt=list(bnd=bnd)
              family="gaussian")
print(summary(F90))
pdf("F90.pdf", useDingbats=FALSE)
  plot(z, fitted.values(F90))
  abline(0,1)
  plot(F90, scheme=2, residuals=TRUE)
  points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
  points(knots$x, knots$y, pch=16, cex=1.1, col=4)
  text(100, 500, "F90")
dev.off()



PGAM <- gamIntervals(predict(F90, se.fit=TRUE), gam=F90, interval="prediction")
#plot(F90$hat, (PGAM$se.fit/PGAM$residual.scale)^2)
#abline(0,1)
sigma <- PGAM$residual.scale[1]
# Terms invert in upper/lower meaning and hence the flipping during data.frame construction.
F90df <- data.frame(comid=Z$comid, site_no=Z$site_no, huc12=Z$huc12,
                     decade=Z$decade,
                     dec_long_va=Z$dec_long_va, dec_lat_va=Z$dec_lat_va,
                     in_model_f90=TRUE, f90=Z$f90,
                     est_lwr_f90=10^PGAM$lwr-1,
                     est_f90    =10^PGAM$fit-1,
                     est_upr_f90=10^PGAM$upr-1, stringsAsFactors=FALSE)
F90df$rse_f90 <- sigma
F90df$se.fit_f90 <- PGAM$se.fit

sites_to_fill <- unique(c(sites_of_area_bust, DDo$site_no[DDo$ed_rch_zone == 1]))
i <- 0
for(site in sites_to_fill) {
  for(decade in as.character(DDo$decade[DDo$site_no == site])) {
    i <- i + 1
    message(site," ", decade, " ", i)
    tmp <- data.frame(comid=DDo$comid[DDo$site_no == site & DDo$decade == decade],
                      site_no=site,
                      huc12=DDo$huc12[DDo$site_no == site & DDo$decade == decade],
                      decade=decade,
                      dec_long_va=DDo$dec_long_va[DDo$site_no == site & DDo$decade == decade],
                      dec_lat_va=DDo$dec_lat_va[DDo$site_no == site & DDo$decade == decade],
                      in_model_f90=FALSE,
                      f90=DDo$f90[DDo$site_no == site & DDo$decade == decade],
                      est_lwr_f90=NA, est_f90=NA, est_upr_f90=NA,
                      rse_f90=NA, se.fit_f90=NA, stringsAsFactors=FALSE)
    F90df <- rbind(F90df, tmp)
  }
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(F90, newdata=tmp, se.fit=TRUE)
}
F90df <- F90df[order(F90df$site_no, F90df$decade),]

for(site in sites_to_fill) {
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(F90, newdata=tmp, se.fit=TRUE)
  pgk <- gamIntervals(jnk, gam=F90, interval="prediction")
  df <- data.frame(comid=tmp$comid, site_no=tmp$site_no, huc12=tmp$huc12,
                   decade=tmp$decade,
                   dec_long_va=tmp$dec_long_va, dec_lat_va=tmp$dec_lat_va,
                   in_model_f90=0,
                   f90=tmp$f90,
                   est_lwr_f90=10^pgk$lwr-1,
                   est_f90    =10^pgk$fit-1,
                   est_upr_f90=10^pgk$upr-1, stringsAsFactors=FALSE)
  df$rse_f90 <- sigma
  df$se.fit_f90 <- pgk$se.fit
  F90df[F90df$site_no == site,] <- df
}
F90df$est_f90[F90df$est_f90 < 0] <- 0
F90df$est_lwr_f90[F90df$est_lwr_f90 < 0] <- 0
F90df$est_upr_f90[F90df$est_upr_f90 < 0] <- 0


#write_feather(F90df, "all_gage_est_f90.feather")



z <- log10(D$f95+1)      # ---------------------------
F95   <- gam(z~s(basin_area) + s(flood_storage, k=7) +
               s(ppt_mean, k=5) + s(dni_ann, k=7)+
               developed+
               mixed_forest+shrubland+
               decade-1+
               s(x, y), knots=knots, data=D, #, bs="so", xt=list(bnd=bnd)
             family="gaussian")
print(summary(F95))
pdf("F95.pdf", useDingbats=FALSE)
plot(z, fitted.values(F95))
abline(0,1)
plot(F95, scheme=2, residuals=TRUE)
points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
points(knots$x, knots$y, pch=16, cex=1.1, col=4)
text(100, 500, "F95")
dev.off()



PGAM <- gamIntervals(predict(F95, se.fit=TRUE), gam=F95, interval="prediction")
#plot(F95$hat, (PGAM$se.fit/PGAM$residual.scale)^2)
#abline(0,1)
sigma <- PGAM$residual.scale[1]
# Terms invert in upper/lower meaning and hence the flipping during data.frame construction.
F95df <- data.frame(comid=Z$comid, site_no=Z$site_no, huc12=Z$huc12,
                     decade=Z$decade,
                     dec_long_va=Z$dec_long_va, dec_lat_va=Z$dec_lat_va,
                     in_model_f95=TRUE, f95=Z$f95,
                     est_lwr_f95=10^PGAM$lwr-1,
                     est_f95    =10^PGAM$fit-1,
                     est_upr_f95=10^PGAM$upr-1, stringsAsFactors=FALSE)
F95df$rse_f95 <- sigma
F95df$se.fit_f95 <- PGAM$se.fit

sites_to_fill <- unique(c(sites_of_area_bust, DDo$site_no[DDo$ed_rch_zone == 1]))
i <- 0
for(site in sites_to_fill) {
  for(decade in as.character(DDo$decade[DDo$site_no == site])) {
    i <- i + 1
    message(site," ", decade, " ", i)
    tmp <- data.frame(comid=DDo$comid[DDo$site_no == site & DDo$decade == decade],
                      site_no=site,
                      huc12=DDo$huc12[DDo$site_no == site & DDo$decade == decade],
                      decade=decade,
                      dec_long_va=DDo$dec_long_va[DDo$site_no == site & DDo$decade == decade],
                      dec_lat_va=DDo$dec_lat_va[DDo$site_no == site & DDo$decade == decade],
                      in_model_f95=FALSE,
                      f95=DDo$f95[DDo$site_no == site & DDo$decade == decade],
                      est_lwr_f95=NA, est_f95=NA, est_upr_f95=NA,
                      rse_f95=NA, se.fit_f95=NA, stringsAsFactors=FALSE)
    F95df <- rbind(F95df, tmp)
  }
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(F95, newdata=tmp, se.fit=TRUE)
}
F95df <- F95df[order(F95df$site_no, F95df$decade),]

for(site in sites_to_fill) {
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(F95, newdata=tmp, se.fit=TRUE)
  pgk <- gamIntervals(jnk, gam=F95, interval="prediction")
  df <- data.frame(comid=tmp$comid, site_no=tmp$site_no, huc12=tmp$huc12,
                   decade=tmp$decade,
                   dec_long_va=tmp$dec_long_va, dec_lat_va=tmp$dec_lat_va,
                   in_model_f95=0,
                   f95=tmp$f95,
                   est_lwr_f95=10^pgk$lwr-1,
                   est_f95    =10^pgk$fit-1,
                   est_upr_f95=10^pgk$upr-1, stringsAsFactors=FALSE)
  df$rse_f95 <- sigma
  df$se.fit_f95 <- pgk$se.fit
  F95df[F95df$site_no == site,] <- df
}
F95df$est_f95[F95df$est_f95 < 0] <- 0
F95df$est_lwr_f95[F95df$est_lwr_f95 < 0] <- 0
F95df$est_upr_f95[F95df$est_upr_f95 < 0] <- 0

#write_feather(F95df, "all_gage_est_f95.feather")




z <- log10(D$f98+1)      # ---------------------------
F98   <- gam(z~s(basin_area) + s(flood_storage, k=7) +
               s(ppt_mean, k=5) + s(dni_ann, k=7)+
               developed+
               mixed_forest+shrubland+
               decade-1+
               s(x, y), knots=knots, data=D, #, bs="so", xt=list(bnd=bnd)
             family="gaussian")
print(summary(F98))
pdf("F98.pdf", useDingbats=FALSE)
plot(z, fitted.values(F98))
abline(0,1)
plot(F98, scheme=2, residuals=TRUE)
points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
points(knots$x, knots$y, pch=16, cex=1.1, col=4)
text(100, 500, "F98")
dev.off()

PGAM <- gamIntervals(predict(F98, se.fit=TRUE), gam=F98, interval="prediction")
#plot(F98$hat, (PGAM$se.fit/PGAM$residual.scale)^2)
#abline(0,1)
sigma <- PGAM$residual.scale[1]
# Terms invert in upper/lower meaning and hence the flipping during data.frame construction.
F98df <- data.frame(comid=Z$comid, site_no=Z$site_no, huc12=Z$huc12,
                     decade=Z$decade,
                     dec_long_va=Z$dec_long_va, dec_lat_va=Z$dec_lat_va,
                     in_model_f98=TRUE, f98=Z$f98,
                     est_lwr_f98=10^PGAM$lwr-1,
                     est_f98    =10^PGAM$fit-1,
                     est_upr_f98=10^PGAM$upr-1, stringsAsFactors=FALSE)
F98df$rse_f98 <- sigma
F98df$se.fit_f98 <- PGAM$se.fit

sites_to_fill <- unique(c(sites_of_area_bust, DDo$site_no[DDo$ed_rch_zone == 1]))
i <- 0
for(site in sites_to_fill) {
  for(decade in as.character(DDo$decade[DDo$site_no == site])) {
    i <- i + 1
    message(site," ", decade, " ", i)
    tmp <- data.frame(comid=DDo$comid[DDo$site_no == site & DDo$decade == decade],
                      site_no=site,
                      huc12=DDo$huc12[DDo$site_no == site & DDo$decade == decade],
                      decade=decade,
                      dec_long_va=DDo$dec_long_va[DDo$site_no == site & DDo$decade == decade],
                      dec_lat_va=DDo$dec_lat_va[DDo$site_no == site & DDo$decade == decade],
                      in_model_f98=FALSE,
                      f98=DDo$f98[DDo$site_no == site & DDo$decade == decade],
                      est_lwr_f98=NA, est_f98=NA, est_upr_f98=NA,
                      rse_f98=NA, se.fit_f98=NA, stringsAsFactors=FALSE)
    F98df <- rbind(F98df, tmp)
  }
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(F98, newdata=tmp, se.fit=TRUE)
}
F98df <- F98df[order(F98df$site_no, F98df$decade),]

for(site in sites_to_fill) {
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(F98, newdata=tmp, se.fit=TRUE)
  pgk <- gamIntervals(jnk, gam=F98, interval="prediction")
  df <- data.frame(comid=tmp$comid, site_no=tmp$site_no, huc12=tmp$huc12,
                   decade=tmp$decade,
                   dec_long_va=tmp$dec_long_va, dec_lat_va=tmp$dec_lat_va,
                   in_model_f98=0,
                   f98=tmp$f98,
                   est_lwr_f98=10^pgk$lwr-1,
                   est_f98    =10^pgk$fit-1,
                   est_upr_f98=10^pgk$upr-1, stringsAsFactors=FALSE)
  df$rse_f98 <- sigma
  df$se.fit_f98 <- pgk$se.fit
  F98df[F98df$site_no == site,] <- df
}
F98df$est_f98[F98df$est_f98 < 0] <- 0
F98df$est_lwr_f98[F98df$est_lwr_f98 < 0] <- 0
F98df$est_upr_f98[F98df$est_upr_f98 < 0] <- 0

#write_feather(F98df, "all_gage_est_f98.feather")




z <- log10(D$f99+1)      # ---------------------------
F99   <- gam(z~s(basin_area) + s(flood_storage, k=7) +
               s(ppt_mean, k=5) + s(dni_ann, k=7)+
               developed+
               mixed_forest+shrubland+
               decade-1+
               s(x, y), knots=knots, data=D, #, bs="so", xt=list(bnd=bnd)
             family="gaussian")
print(summary(F99))
pdf("F99.pdf", useDingbats=FALSE)
plot(z, fitted.values(F99))
abline(0,1)
plot(F99, scheme=2, residuals=TRUE)
points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
points(knots$x, knots$y, pch=16, cex=1.1, col=4)
text(100, 500, "F99")
dev.off()


PGAM <- gamIntervals(predict(F99, se.fit=TRUE), gam=F99, interval="prediction")
#plot(F99$hat, (PGAM$se.fit/PGAM$residual.scale)^2)
#abline(0,1)
sigma <- PGAM$residual.scale[1]
# Terms invert in upper/lower meaning and hence the flipping during data.frame construction.
F99df <- data.frame(comid=Z$comid, site_no=Z$site_no, huc12=Z$huc12,
                     decade=Z$decade,
                     dec_long_va=Z$dec_long_va, dec_lat_va=Z$dec_lat_va,
                     in_model_f99=TRUE, f99=Z$f99,
                     est_lwr_f99=10^PGAM$lwr-1,
                     est_f99    =10^PGAM$fit-1,
                     est_upr_f99=10^PGAM$upr-1, stringsAsFactors=FALSE)
F99df$rse_f99 <- sigma
F99df$se.fit_f99 <- PGAM$se.fit

sites_to_fill <- unique(c(sites_of_area_bust, DDo$site_no[DDo$ed_rch_zone == 1]))
i <- 0
for(site in sites_to_fill) {
  for(decade in as.character(DDo$decade[DDo$site_no == site])) {
    i <- i + 1
    message(site," ", decade, " ", i)
    tmp <- data.frame(comid=DDo$comid[DDo$site_no == site & DDo$decade == decade],
                      site_no=site,
                      huc12=DDo$huc12[DDo$site_no == site & DDo$decade == decade],
                      decade=decade,
                      dec_long_va=DDo$dec_long_va[DDo$site_no == site & DDo$decade == decade],
                      dec_lat_va=DDo$dec_lat_va[DDo$site_no == site & DDo$decade == decade],
                      in_model_f99=FALSE,
                      f99=DDo$f99[DDo$site_no == site & DDo$decade == decade],
                      est_lwr_f99=NA, est_f99=NA, est_upr_f99=NA,
                      rse_f99=NA, se.fit_f99=NA, stringsAsFactors=FALSE)
    F99df <- rbind(F99df, tmp)
  }
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(F99, newdata=tmp, se.fit=TRUE)
}
F99df <- F99df[order(F99df$site_no, F99df$decade),]

for(site in sites_to_fill) {
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(F99, newdata=tmp, se.fit=TRUE)
  pgk <- gamIntervals(jnk, gam=F99, interval="prediction")
  df <- data.frame(comid=tmp$comid, site_no=tmp$site_no, huc12=tmp$huc12,
                   decade=tmp$decade,
                   dec_long_va=tmp$dec_long_va, dec_lat_va=tmp$dec_lat_va,
                   in_model_f99=0,
                   f99=tmp$f99,
                   est_lwr_f99=10^pgk$lwr-1,
                   est_f99    =10^pgk$fit-1,
                   est_upr_f99=10^pgk$upr-1, stringsAsFactors=FALSE)
  df$rse_f99 <- sigma
  df$se.fit_f99 <- pgk$se.fit
  F99df[F99df$site_no == site,] <- df
}
F99df$est_f99[F99df$est_f99 < 0] <- 0
F99df$est_lwr_f99[F99df$est_lwr_f99 < 0] <- 0
F99df$est_upr_f99[F99df$est_upr_f99 < 0] <- 0

#write_feather(F99df, "all_gage_est_f99.feather")



z <- log10(D$f99.9+1)      # ---------------------------
F99p9   <- gam(z~s(basin_area) + s(flood_storage, k=7) +
               s(ppt_mean, k=5) + s(dni_ann, k=7)+
               developed+
               mixed_forest+shrubland+
               decade-1+
               s(x, y), knots=knots, data=D, #, bs="so", xt=list(bnd=bnd)
             family="gaussian")
print(summary(F99p9))
pdf("F99p9.pdf", useDingbats=FALSE)
plot(z, fitted.values(F99p9))
abline(0,1)
plot(F99p9, scheme=2, residuals=TRUE)
points(D$x, D$y, pch=4, lwd=.5, cex=0.9, col=8)
points(knots$x, knots$y, pch=16, cex=1.1, col=4)
text(100, 500, "F99p9")
dev.off()



PGAM <- gamIntervals(predict(F99p9, se.fit=TRUE), gam=F99p9, interval="prediction")
#plot(F99p9$hat, (PGAM$se.fit/PGAM$residual.scale)^2)
#abline(0,1)
sigma <- PGAM$residual.scale[1]
# Terms invert in upper/lower meaning and hence the flipping during data.frame construction.
F99p9df <- data.frame(comid=Z$comid, site_no=Z$site_no, huc12=Z$huc12,
                     decade=Z$decade,
                     dec_long_va=Z$dec_long_va, dec_lat_va=Z$dec_lat_va,
                     in_model_f99.9=TRUE, f99.9=Z$f99.9,
                     est_lwr_f99.9=10^PGAM$lwr-1,
                     est_f99.9    =10^PGAM$fit-1,
                     est_upr_f99.9=10^PGAM$upr-1, stringsAsFactors=FALSE)
F99p9df$rse_f99.9 <- sigma
F99p9df$se.fit_f99.9 <- PGAM$se.fit

sites_to_fill <- unique(c(sites_of_area_bust, DDo$site_no[DDo$ed_rch_zone == 1]))
i <- 0
for(site in sites_to_fill) {
  for(decade in as.character(DDo$decade[DDo$site_no == site])) {
    i <- i + 1
    message(site," ", decade, " ", i)
    tmp <- data.frame(comid=DDo$comid[DDo$site_no == site & DDo$decade == decade],
                      site_no=site,
                      huc12=DDo$huc12[DDo$site_no == site & DDo$decade == decade],
                      decade=decade,
                      dec_long_va=DDo$dec_long_va[DDo$site_no == site & DDo$decade == decade],
                      dec_lat_va=DDo$dec_lat_va[DDo$site_no == site & DDo$decade == decade],
                      in_model_f99.9=FALSE,
                      f99.9=DDo$f99.9[DDo$site_no == site & DDo$decade == decade],
                      est_lwr_f99.9=NA, est_f99.9=NA, est_upr_f99.9=NA,
                      rse_f99.9=NA, se.fit_f99.9=NA, stringsAsFactors=FALSE)
    F99p9df <- rbind(F99p9df, tmp)
  }
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(F99p9, newdata=tmp, se.fit=TRUE)
}
F99p9df <- F99p9df[order(F99p9df$site_no, F99p9df$decade),]

for(site in sites_to_fill) {
  tmp <- DDo[DDo$site_no == site,]
  jnk <- predict(F99p9, newdata=tmp, se.fit=TRUE)
  pgk <- gamIntervals(jnk, gam=F99p9, interval="prediction")
  df <- data.frame(comid=tmp$comid, site_no=tmp$site_no, huc12=tmp$huc12,
                   decade=tmp$decade,
                   dec_long_va=tmp$dec_long_va, dec_lat_va=tmp$dec_lat_va,
                   in_model_f99.9=0,
                   f99.9=tmp$f99.9,
                   est_lwr_f99.9=10^pgk$lwr-1,
                   est_f99.9    =10^pgk$fit-1,
                   est_upr_f99.9=10^pgk$upr-1, stringsAsFactors=FALSE)
  df$rse_f99.9 <- sigma
  df$se.fit_f99.9 <- pgk$se.fit
  F99p9df[F99p9df$site_no == site,] <- df
}
F99p9df$est_f99.9[F99p9df$est_f99.9 < 0] <- 0
F99p9df$est_lwr_f99.9[F99p9df$est_lwr_f99.9 < 0] <- 0
F99p9df$est_upr_f99.9[F99p9df$est_upr_f99.9 < 0] <- 0

#write_feather(F99p9df, "all_gage_est_f99p9.feather")

#sink()

save(F50, F90, F95, F98, F99, F99p9, file="QuantileModels.RData")



plot(log10(F50df$f50), log10(F50df$est_f50),
     xlab="Observed value", ylab="Fitted value")
mtext(paste0("f50: NSE=",round(NSE(F50df$est_f50, F50df$f50), digits=2)))
abline(0,1)

plot(log10(F90df$f90), log10(F90df$est_f90),
     xlab="Observed value", ylab="Fitted value")
mtext(paste0("f90: NSE=",round(NSE(F90df$est_f90, F90df$f90), digits=2)))
abline(0,1)

plot(log10(F95df$f95), log10(F95df$est_f95),
     xlab="Observed value", ylab="Fitted value")
mtext(paste0("f95: NSE=",round(NSE(F95df$est_f95, F95df$f95), digits=2)))
abline(0,1)

plot(log10(F98df$f98), log10(F98df$est_f98),
     xlab="Observed value", ylab="Fitted value")
mtext(paste0("f98: NSE=",round(NSE(F98df$est_f98, F98df$f98), digits=2)))
abline(0,1)

plot(log10(F99df$f99), log10(F99df$est_f99),
     xlab="Observed value", ylab="Fitted value")
mtext(paste0("f99: NSE=",round(NSE(F99df$est_f99, F99df$f99), digits=2)))
abline(0,1)

plot(log10(F99p9df$f99.9), log10(F99p9df$est_f99.9),
     xlab="Observed value", ylab="Fitted value")
mtext(paste0("f99.9: NSE=",round(NSE(F99p9df$est_f99.9, F99p9df$f99.9), digits=2)))
abline(0,1)

