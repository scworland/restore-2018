library(akqdecay); library(sp); library(lmomco)

load("DV.RData")

TFDC <- new.env()
n <- fill_tfdcenv(dvenv=DV, envir=TFDC) # number fdctrend() processed
AllFDCtrend <- visFDCtrend(TFDC, file="fdc_trend.pdf")
save(AllFDCtrend, TFDC, file="AllFDCtrend.RData")
#visFDCtrend(TFDC$"08093500") # new feature 09/28/2017

# 12/01/2017
sitesWOcompleteFDC <- c("02295163", "02300082", "02388350", "02432500",
                        "07031692", "08020450", "08067000", "08106350")


FDCLMRlog <- new.env()
fill_lmrfdcenv(dvenv=DV, envir=FDCLMRlog, log=TRUE)
FDCLMRnolog <- new.env()
fill_lmrfdcenv(dvenv=DV, envir=FDCLMRnolog, log=FALSE)

FDClmrdf.log   <- lmrfdc_table(FDCLMRlog)
FDClmrdf.nolog <- lmrfdc_table(FDCLMRnolog)

save(FDCLMRlog, FDCLMRnolog, FDClmrdf.log, FDClmrdf.nolog, file="FDCLMR.RData")
load("FDCLMR.RData")

pdf("lmrdia.pdf", useDingbats=FALSE)
  L <- FDClmrdf.log
  L$color <- rgb(.5,.3,0,.3)
  L$color[L$nzero > 30 ] <- rgb(.3, 0.8, 0.1, 0.5) # greenish towards
  L$color[L$nzero > 60 ] <- rgb(.3, 0.6, 0.2, 0.5)
  L$color[L$nzero > 90 ] <- rgb(.3, 0.4, 0.3, 0.5)
  L$color[L$nzero > 120] <- rgb(.3, 0.2, 0.4, 0.5)
  L$color[L$nzero > 180] <- rgb(.3, 0,.5,0.5) # purple

  plotlmrdia(lmrdia(), xlim=c(-.6,.9), ylim=c(-.2,0.6))
  points(L$T3, L$T4, cex=0.3, col=L$color, pch=16, lwd=0.3)
  mtext("log-transform of (DV+1cfs)")

  N <- FDClmrdf.nolog
  plotlmrdia(lmrdia(), xlim=c(-.4,1), ylim=c(-.2,1))
  points(N$T3, N$T4, cex=0.3, col=rgb(.5,.3,0,.3), pch=16, lwd=0.3)
  mtext("no transformation")
dev.off()


stop()  # experimental material follows below

for(ktsites in unique(AllFDCtrend$site)) {
   X <- AllFDCtrend[AllFDCtrend$site == ktsites,]
   message(paste0(ktsites," ",length(X$site)))
   deltau <- 0.05; alpha <- 0.05
   rnames <- seq(-1,1,by=deltau)
   cnames <- round((1:365)/(365+1), digits=3)
   nr <- length(rnames); nc <- length(cnames)
   MX <- CX <- TX <- matrix(NA, nrow=nr, ncol=nc)
   rownames(MX) <- rnames; rownames(CX) <- rnames
   colnames(MX) <- cnames; colnames(CX) <- cnames
   rix <- 1:length(rnames)
   for(column in 1:365) {
      tau <- X$estimate[column]
      color <- as.numeric(X$p.value[column] > alpha); color[color != 1] <- 2
      ix <- sample(rix[rnames >= tau-deltau & rnames <= tau+deltau], 1)
      TX[1:ix, column] <- 3; TX[ix:nr,column] <- 4
      MX[ix,column] <- tau
      CX[ix,column] <- TX[ix,column] <- color
   }
  stop()
}
#jpeg(filename="junkMX.jpg", width=365, height=401)
#image(t(MX), col=topo.colors(1), xaxt="n", yaxt="n")
#dev.off()
#jpeg(filename="junkCX.jpg", width=365, height=401)
#image(t(CX), col=topo.colors(2), xaxt="n", yaxt="n")
#dev.off()
