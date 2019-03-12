library(sp)
library(rgdal)

LATLONG <- paste0("+proj=longlat +ellps=GRS80 ",
                  "+datum=NAD83 +no_defs +towgs84=0,0,0")
LATLONG <- CRS(LATLONG)
ALBEA <- paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 ",
                "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
ALBEA <- CRS(ALBEA)

SB   <- read.table("../scibase/huc12/all_huc12_covariates.csv", sep=",", header=TRUE,
                  colClasses="character")
SB   <- SpatialPointsDataFrame(cbind(as.numeric(SB$dec_long_va),
                                     as.numeric(SB$dec_lat_va)),
                                data=SB, proj4string=LATLONG)
SB <- spTransform(SB, ALBEA)
SB <- SB[SB$decade == "1950",]

load("./restore_huc12.RData")
F12 <- huc12

F12@data$HUC12 <- as.character(F12@data$HUC12)
R12 <- data.frame(HUC12=unique(SB$huc12), stringsAsFactors=FALSE)
A12 <- data.frame(HUC12=F12@data$HUC12,   stringsAsFactors=FALSE)
F12@data$IN_COVARIATES <- FALSE

not.found <- NA
for(huc in unique(SB$huc12)) {
  if(length(F12@data$HUC12[F12@data$HUC12 == huc]) == 1) {
    #message(huc," is found")
    F12@data$IN_COVARIATES[F12@data$HUC12 == huc] <- TRUE
  } else {
    message(huc," is not found")
    not.found <- c(not.found, huc)
  }
}
not.found <- sort(not.found[! is.na(not.found)])
SB$ECO_SB <- FALSE

pdf("testHUC12.pdf", useDingbats=FALSE, width=11, height=10)
  plot(F12[  F12@data$IN_COVARIATES,], col=4, lwd=0.1)
  plot(F12[! F12@data$IN_COVARIATES,], col=2, lwd=0.1, add=TRUE)
  for(huc in not.found) {
    #message("plotting ", huc)
    #tmp <- SB[SB$huc12 == huc,]
    #plot(tmp, col=2, add=TRUE, pch=16, cex=0.4)
    SB$ECO_SB[SB$huc12 == huc] <- 1
  }
  JK <- SB[SB$ECO_SB == TRUE,]
  plot(JK, add=TRUE, pch=21, col=6, bg=3)
dev.off()


