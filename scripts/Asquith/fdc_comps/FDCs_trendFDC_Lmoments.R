library(akqdecay); library(lmomco)

load("DV.RData")

TFDC <- new.env()
n <- fill_tfdcenv(dvenv=DV, envir=TFDC) # number fdctrend() processed
AllFDCtrend <- visFDCtrend(TFDC, file="all_fdc_trend.pdf")
save(AllFDCtrend, TFDC, file="AllFDCtrend.RData")
#visFDCtrend(TFDC$"08093500", site="08093500") # new feature 09/28/2017

# 02/15/2018
sitesWOcompleteFDC <- c("02295163", "02388350", "02432500",
                        "07031692", "08020450", "08067000", "08106350")


FDCLMRlog <- new.env()
fill_lmrfdcenv(dvenv=DV, envir=FDCLMRlog,   log=TRUE,  decade=TRUE)
FDCLMRnolog <- new.env()
fill_lmrfdcenv(dvenv=DV, envir=FDCLMRnolog, log=FALSE, decade=TRUE)

FDClmrdf.log   <- lmrfdc_table(FDCLMRlog)
FDClmrdf.nolog <- lmrfdc_table(FDCLMRnolog)

save(FDCLMRlog, FDCLMRnolog, FDClmrdf.log, FDClmrdf.nolog, file="FDCLMR.RData")
load("FDCLMR.RData")

pdf("lmrdia.pdf", useDingbats=FALSE)
  L <- FDClmrdf.log; N <- FDClmrdf.nolog
  L$color <- N$color <- rgb(.5,.3,0,.3) # zeros will go greenish towards purple
  L$color[L$nzero > 30 ] <- N$color[N$nzero > 30 ] <- rgb(.3,.8,.1,.5) # green
  L$color[L$nzero > 60 ] <- N$color[N$nzero > 60 ] <- rgb(.3,.6,.2,.5)
  L$color[L$nzero > 90 ] <- N$color[N$nzero > 90 ] <- rgb(.3,.4,.3,.5)
  L$color[L$nzero > 120] <- N$color[N$nzero > 120] <- rgb(.3,.2,.4,.5)
  L$color[L$nzero > 180] <- N$color[N$nzero > 180] <- rgb(.3,.0,.5,.5) # purple

  plotlmrdia(lmrdia(), xlim=c(-.5,.9), ylim=c(-.1,0.6))
  points(L$T3, L$T4, cex=0.3, col=L$color, pch=16, lwd=0.3)
  mtext("log-transform of (DV+1cfs)")

  plotlmrdia(lmrdia(), xlim=c(-.05,1), ylim=c(0,1))
  points(N$T3, N$T4, cex=0.3, col=N$color, pch=16, lwd=0.3)
  mtext("no transformation")
dev.off()

modFDClmrdf.log <- FDClmrdf.log[! is.na(FDClmrdf.log$site), ]
modFDClmrdf.log <- modFDClmrdf.log[! is.na(modFDClmrdf.log$min),  ]
modFDClmrdf.log <- modFDClmrdf.log[modFDClmrdf.log$decade >= 1950,]

modFDClmrdf.nolog <- FDClmrdf.nolog[! is.na(FDClmrdf.nolog$site),]
modFDClmrdf.nolog <- modFDClmrdf.nolog[! is.na(modFDClmrdf.nolog$min), ]
modFDClmrdf.nolog <- modFDClmrdf.nolog[modFDClmrdf.nolog$decade >= 1950,]

ModelSitesA <- unique(modFDClmrdf.log$site)
ModelSitesB <- unique(modFDClmrdf.nolog$site)
if(length(ModelSitesA) != length(ModelSitesB)) {
   stop("SOMETHING IS WRONG, STOP THE PROJECT!")
}

ModelSites <- ModelSitesA; rm(ModelSitesA); rm(ModelSitesB)

ModTFDC <- list2env(akq_rm(ModelSites, TFDC, invert=TRUE))
ModFDCtrend <- visFDCtrend(ModTFDC, file="mod_fdc_trend.pdf")
save(ModTFDC, ModFDCtrend, file="ModFDCtrend.RData")


modFDCLMRlog   <- list2env(akq_rm(ModelSites, FDCLMRlog,   invert=TRUE))
modFDCLMRnolog <- list2env(akq_rm(ModelSites, FDCLMRnolog, invert=TRUE))

if(length(ls(modFDCLMRlog)) != length(ls(modFDCLMRnolog)) ) {
   stop("SOMETHING IS WRONG, STOP THE PROJECT!")
}

save(modFDCLMRlog, modFDCLMRnolog, modFDClmrdf.log, modFDClmrdf.nolog, file="modFDCLMR.RData")
load("modFDCLMR.RData")


library(feather)
fdc_lmr_pplo <- modFDClmrdf.log
names <- names(fdc_lmr_pplo)
names[1] <- "site_no"; names(fdc_lmr_pplo) <- names
write_feather(fdc_lmr_pplo, "log_fdc_lmr_pplo.feather")
fdc_lmr_pplo <- read_feather("log_fdc_lmr_pplo.feather")


fdc_lmr_pplo <- modFDClmrdf.nolog
names <- names(fdc_lmr_pplo)
names[1] <- "site_no"; names(fdc_lmr_pplo) <- names
write_feather(fdc_lmr_pplo, "fdc_lmr_pplo.feather")
fdc_lmr_pplo <- read_feather("fdc_lmr_pplo.feather")



ModelSites <- data.frame(site_no=ModelSites)
write.table(ModelSites, file="decade1950plus_site_list.csv", row.names=FALSE)

pdf("mod_lmrdia.pdf", useDingbats=FALSE)
  L <- modFDClmrdf.log; N <- modFDClmrdf.nolog
  L$color <- N$color <- rgb(.5,.3,0,.3) # zeros will go greenish towards purple
  L$color[L$nzero > 30 ] <- N$color[N$nzero > 30 ] <- rgb(.3,.8,.1,.5) # green
  L$color[L$nzero > 60 ] <- N$color[N$nzero > 60 ] <- rgb(.3,.6,.2,.5)
  L$color[L$nzero > 90 ] <- N$color[N$nzero > 90 ] <- rgb(.3,.4,.3,.5)
  L$color[L$nzero > 120] <- N$color[N$nzero > 120] <- rgb(.3,.2,.4,.5)
  L$color[L$nzero > 180] <- N$color[N$nzero > 180] <- rgb(.3,.0,.5,.5) # purple

  plotlmrdia(lmrdia(), xlim=c(-.5,.9), ylim=c(-.1,0.6))
  points(L$T3, L$T4, cex=0.3, col=L$color, pch=16, lwd=0.3)
  mtext("log-transform of (DV+1cfs)")

  plotlmrdia(lmrdia(), xlim=c(-.05,1), ylim=c(0,1))
  points(N$T3, N$T4, cex=0.3, col=N$color, pch=16, lwd=0.3)
  mtext("no transformation")
dev.off()



stop()  # experimental material follows below
