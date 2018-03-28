library(lmomco)

"getSB" <- function(sp.points, ...) {
  message("Select points and FINISH (RStudio)")
  xy <- coordinates(sp.points)
  dat <- slot(sp.points, "data")
  whos <- identify(xy[,1], xy[,2], cex=0.8, col=grey(.1), ...)
  for(who in whos) {
     message("site=",who," is '",   dat$SITE_BADGE[who])
  }
  return(whos)
}


spCOV$aep4x <- NA
spCOV$aep4a <- NA
spCOV$aep4k <- NA
spCOV$aep4h <- NA
n <- length(spCOV$comid); fails <- rep(NA, n)
for(i in 1:n) {
  message(i,"(",n,")")
  tmp <- spCOV[i,]; rats <- c(NA, tmp$est_T2, tmp$est_T3, tmp$est_T4)
  t34 <- (5*rats[3]^2)/4
  if(rats[3] >  1)  rats[3] <-  0.99
  if(rats[4] >  1)  rats[4] <-  0.99
  if(rats[3] < -1)  rats[3] <- -0.99
  if(rats[4] < -1)  rats[4] <- -0.99
  if(rats[4] < t34) rats[4] <-   t34
  suppressWarnings(lmr <- vec2lmom(c(10^tmp$est_L1,rats[2:4]), lscale=FALSE))
  aep4 <- paraep4(lmr, snap.tau4=TRUE, checklmom = FALSE)
  ifail <- aep4$ifail
  if(ifail) fails[i] <- i
  spCOV$aep4x[i] <- aep4$para[1]
  spCOV$aep4a[i] <- aep4$para[2]
  spCOV$aep4k[i] <- aep4$para[3]
  spCOV$aep4h[i] <- aep4$para[4]
  if(is.na(spCOV$aep4x[i])) stop("failure in aep4")
}

save(spCOV, file="spCOV.RData")




spCOV$gnox  <- NA
spCOV$gnoa  <- NA
spCOV$gnok  <- NA
n <- length(spCOV$comid); fails <- rep(NA, n)
for(i in 1:n) {
  message(i,"(",n,")")
  tmp <- spCOV[i,]; rats <- c(NA, tmp$est_T2, tmp$est_T3, tmp$est_T4)
  if(rats[3] >  0.95)  rats[3] <-  0.95
  if(rats[4] >  0.95)  rats[4] <-  0.95
  if(rats[3] < -0.95)  rats[3] <- -0.95
  if(rats[4] < -0.95)  rats[4] <- -0.95
  if(rats[4] < t34) rats[4] <-   t34
  suppressWarnings(lmr <- vec2lmom(c(L1$duan_smearing*10^tmp$est_L1,rats[2:4]), lscale=FALSE))
  gno <- pargno(lmr, checklmom = FALSE)
  spCOV$gnox[i] <- gno$para[1]
  spCOV$gnoa[i] <- gno$para[2]
  spCOV$gnok[i] <- gno$para[3]
}

save(spCOV, file="spCOV.RData")



spCOV$kapx  <- NA
spCOV$kapa  <- NA
spCOV$kapk  <- NA
spCOV$kaph  <- NA
n <- length(spCOV$comid); fails <- rep(NA, n)
for(i in 1:n) {
  message(i,"(",n,")")
  tmp <- spCOV[i,]; rats <- c(NA, tmp$est_T2, tmp$est_T3, tmp$est_T4)
  t34 <- (5*rats[3]^2)/4
  if(rats[3] >  1)  rats[3] <-  0.99
  if(rats[4] >  1)  rats[4] <-  0.99
  if(rats[3] < -1)  rats[3] <- -0.99
  if(rats[4] < -1)  rats[4] <- -0.99
  if(rats[4] < t34) rats[4] <-   t34
  suppressWarnings(lmr <- vec2lmom(c(10^tmp$est_L1,rats[2:4]), lscale=FALSE))
  kap <- parkap2(lmr, snap.tau4=TRUE, checklmom = FALSE)
  ifail <- kap$ifail
  if(ifail) fails[i] <- i
  spCOV$kapx[i] <- kap$para[1]
  spCOV$kapa[i] <- kap$para[2]
  spCOV$kapk[i] <- kap$para[3]
  spCOV$kaph[i] <- kap$para[4]
  if(is.na(spCOV$kapx[i])) stop("failure in kap")
}

save(spCOV, file="spCOV.RData")





spCOV$gldx  <- NA
spCOV$glda  <- NA
spCOV$gldk  <- NA
spCOV$gldh  <- NA
n <- length(spCOV$comid)
for(i in 1:10) {
  message(i,"(",n,")")
  tmp <- spCOV[i,]; rats <- c(NA, tmp$est_T2, tmp$est_T3, tmp$est_T4, tmp$est_T5)
  t34 <- (5*rats[3]^2)/4
  if(rats[3] >  1)  rats[3] <-  0.99
  if(rats[4] >  1)  rats[4] <-  0.99
  if(rats[5] >  1)  rats[5] <-  0.99
  if(rats[3] < -1)  rats[3] <- -0.99
  if(rats[4] < -1)  rats[4] <- -0.99
  if(rats[5] < -1)  rats[5] <- -0.99
  if(rats[4] < t34) rats[4] <-   t34
  suppressWarnings(lmr <- vec2lmom(c(tmp$est_L1,rats[2:5]), lscale=FALSE))
  gld <- pargld(lmr, checklmom = FALSE)
  next
  spCOV$gldx[i] <- gld$para[1]
  spCOV$glda[i] <- gld$para[2]
  spCOV$gldk[i] <- gld$para[3]
  spCOV$gldh[i] <- gld$para[4]
  if(is.na(spCOV$gldx[i])) stop("failure in gld")
}




load("spCOV.Rdata")
load("DEMO.RData")
load("Models.RData")
library(lmomco)
pdf("test_fdc.pdf", useDingbats=FALSE, height=6, width=7)
opts <- par(); par(las=1, mgp=c(3,0.5,0))
FF <- 1/(3653+1); FF <- seq(qnorm(FF), qnorm(1-FF), 0.005); FF <- pnorm(FF)
comid <- sample(spCOV$comid, 1)
comid <- "3588922" # at Guad Comfort
COV <- spCOV@data
cmp <- COV[COV$comid == comid,]
cmp_area <- 10^cmp$acc_basin_area[1]
plot(1,1, xlim=range(qnorm(FF)), ylim=c(0.1,10000), xaxt="n", yaxt="n", xlab="",
          log="y", type="n",
          ylab="COMID-scaled, daily-mean streamflow, cubic feet per second")
add.lmomco.axis(las=2, tcl=0.5, side.type="NPP", cex=0.8, case="lower")
add.log.axis(side=2, tcl=0.8*abs(par()$tcl), two.sided=TRUE)
add.log.axis(logs=c(1),   tcl=-0.5*abs(par()$tcl), side=2, two.sided=TRUE)
add.log.axis(logs=c(1),   tcl=+1.3*abs(par()$tcl), side=2, two.sided=TRUE)
add.log.axis(logs=c(1, 2, 4, 6), side=2, make.labs=TRUE, las=1, label="")



gage <- D[D$site_no == "08167000",]
gage_area <- 10^gage$acc_basin_area[1]
fdc <- as.data.frame(gage[,c(11:37)])
fdc$coords.x1 <- NULL
fdc$coords.x2 <- NULL
ff <- as.numeric(gsub("f","",names(fdc)))/100
ff <- qnorm(ff)
fdc <- as.data.frame(t(as.matrix(fdc)))
names(fdc) <- as.character(gage$decade)
phi <- 0.9 # Asquith and others (2006): Drainage-area ratio
col <- .8; pch <- 3
for(i in 1:length(names(fdc))) {
  gmp <- fdc[,i]; gmp[gmp == 0] <- 0.001
  lines(ff, gmp*(cmp_area/gage_area)^phi, col=grey(col))
  points(ff, gmp*(cmp_area/gage_area)^phi, pch=pch, col=grey(col))
  col <- col - .1; pch <- ifelse(pch == 3, 4, 3)
}
text(-3.5, 10000, "USGS 08167000 streamgage ('next to' COMID 3588922) has area of 2,181 square kilometers", pos=4, cex=0.7)
text(-3.5, 7000, "COMID 3588922 has area of 189 square kilometers", pos=4, cex=0.7)
text(-3.5, 3000, "'Observed' flow-duration curve (FDC) for COMID = \n Gaged FDC * (acc_basin_area of COMID / acc_basin_area of streamgage)^0.9", pos=4, cex=0.7)
text(-1.0,.15, paste0("A Duan smearing estimator (Helsel and Hirsch, 2002) of ",
                     round(L1$duan_smearing, digits=3), " was\n",
                     "used to correct retransformation bias in the mean (first L-moment)\n",
                     "prior to fitting of the distributions to the L-moments."),
             pos=4, cex=0.7)
txt <- c(paste0("Observed flow-duration curve by decade where grey\n",
                "increases from 1950s to 2000s and\nsymbol alternates by decade"))
legend(-3.5, 2000, txt, col=grey(0.2), bty="n", lwd=1, pch=NA, cex=0.7)



for(decade in cmp$decade) {
  lwd <- ifelse(decade == "1950", 0.5,
         ifelse(decade == "1960", 1,
         ifelse(decade == "1970", 1.5,
         ifelse(decade == "1980", 0.5,
         ifelse(decade == "1990", 1.0,
         ifelse(decade == "2000", 1.5, 2))))))
  lty <- ifelse(decade == "1950", 2,
         ifelse(decade == "1960", 2,
         ifelse(decade == "1970", 2,
         ifelse(decade == "1980", 1,
         ifelse(decade == "1990", 1,
         ifelse(decade == "2000", 1, 3))))))
  message(decade)
  tmp <- cmp[cmp$decade == decade,]
  tmp <- tmp[1,]
  aep4 <- c(tmp$aep4x,tmp$aep4a, tmp$aep4k,tmp$aep4h)
  aep4 <- vec2par(aep4, type="aep4")
  nep <- qnorm(f2f(FF, pp=tmp$est_pplo))
  qua <- qlmomco(f2flo(FF, pp=tmp$est_pplo), aep4)
  qua[qua <= 0.01] <- 0.01
  lines(nep,qua, col=rgb(1,0,0,0.9), lwd=lwd, lty=lty)
  gno <- c(tmp$gnox,tmp$gnoa, tmp$gnok)
  gno <- vec2par(gno, type="gno")
  nep <- qnorm(f2f(FF, pp=tmp$est_pplo))
  qua <- qlmomco(f2flo(FF, pp=tmp$est_pplo), gno)
  qua[qua <= 0.01] <- 0.01
  lines(nep,qua, col=rgb(0,1,0,0.9), lwd=lwd, lty=lty)
  kap <- c(tmp$kapx,tmp$kapa, tmp$kapk,tmp$kaph)
  kap <- vec2par(kap, type="kap")
  nep <- qnorm(f2f(FF, pp=tmp$est_pplo))
  qua <- qlmomco(f2flo(FF, pp=tmp$est_pplo), kap)
  qua[qua <= 0.01] <- 0.01
  lines(nep,qua, col=rgb(0,0,1,0.9), lwd=lwd, lty=lty)
}
mtext(paste0("COMID: ",comid))
txt <- c("Asymmetric exponential power distribution (four parameter)",
         "Generalized normal distribution (three parameter log-Normal)",
         "Kappa distribution (four parameter)")
legend(-3.5, 300, txt, col=c(2,3,4), bty="n", lwd=1, cex=0.7)
suppressWarnings(par(opts))
dev.off()

FF <- c(0,FF,1)
new_qua <- c(rep(0, length(FF)-length(qua)), qua)
