library(sp)
library(rgeos)
library(rgdal)
library(GISTools)
library(feather)

load(file.choose()) # RESTOREstreams.RData
load(file.choose()) # GulfStates.RData
load(file.choose()) # spRESTORE_MGCV_BND.RData
load(file.choose()) # Edwards

DD <- read_feather(file.choose())

LATLONG <- paste0("+init=epsg:4269 +proj=longlat +ellps=GRS80 ",
                  "+datum=NAD83 +no_defs +towgs84=0,0,0")
LATLONG <- sp::CRS(LATLONG)
ALBEA <- paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 ",
                "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
ALBEA <- sp::CRS(ALBEA)

DD <- SpatialPointsDataFrame(cbind(DD$lon, DD$lat), data=DD, proj4string=LATLONG)
DD <- spTransform(DD, ALBEA)


east_grids  <- seq(80,100,by=2)
north_grids <- seq(26,38, by=2)

gx <- gy <- vector(mode="numeric")
for(i in 1:length(north_grids)) {
   gy <- c(gy, rep(north_grids[i], length(east_grids)))
   gx <- c(gx, east_grids)
}
GL <- SpatialPoints(cbind(-gx, gy), proj4string=LATLONG)
GL <- spTransform(GL, ALBEA)
XY <- coordinates(GL)
x <- XY[,1]; y <- XY[,2]
#ind <- mgcv::inSide(bnd,x,y)
#XY <- XY[ind,]
GL <- SpatialPoints(cbind(x,y), proj4string=ALBEA)
ix <- 1:length(x)
plot(GL, pch=1, col=2)
text(XY[,1],XY[,2], ix)
GL <- GL[-c(1,3:9, 12, 14:20, 23, 34, 45)]



map_annotation <- function() {
  ss <- list(x=-90000, y=320000)
  SpatialPolygonsRescale(layout.scale.bar(height=.07), offset=ss, scale=200*1000,
                         fill=c("transparent", "black"), plot.grid=FALSE, lwd=0.4)
  ss <- list(x=-90000, y=303000)
  SpatialPolygonsRescale(layout.scale.bar(height=.08), offset=ss, scale=100*1609.344,
                         fill=c("transparent", "black"), plot.grid=FALSE, lwd=0.4)
  xx <- -122000
  sl <- list(x=xx, y=350000)
  sr <- list(x=xx+200*1000, y=350000)
  text(sl, "0",xx, cex=0.6, pos=4); text(sr, "200 kilometers", cex=0.5, pos=4)
  sl <- list(x=xx, y=283000)
  sr <- list(x=xx+100*1609.344, y=283000)
  text(sl, "0", cex=0.6, pos=4); text(sr, "100 miles", cex=0.5, pos=4)

  txt <- paste0("Albers Equal Area Projection\n",
                "North American Datum of 1983\n",
                "Base modified from U.S. Geological Survey digital data, 1:24,000")
  text(-420000, 1590000, txt, pos=4, cex=0.45)
  #plot(GulfStates_modified, add=TRUE, lwd=.4, lty=2)
  STATES <- c("Texas", "Oklahoma", "Missouri", "Arkansas", "Louisiana", "Mississippi",
              "Tennessee", "Kentucky", "Alabama", "Georgia", "Florida")
  STATES <- data.frame(easting=c(-440000, -202900.4,  178000,  178000,  279961.5,
                                 430000,  710000,  710000,  710000,
                                 1100000, 1270000),
                       northing=c(955139.0, 1400000, 1558716.4, 1400000,  795368.5,
                                  1170000, 1450000, 1600000, 1325000,
                                  1170000, 800000),
                       state=STATES)
  text(STATES$easting, STATES$northing, STATES$state, pos=4, cex=0.7, col=grey(0.3))
  plot(GL, lwd=0.4, col=grey(0.22), add=TRUE)
}


pdf("SiteMap_junk.pdf", useDingbats=FALSE, width=11, height=10)
  plot(spRESTORE_MGCV_BND)
  usr <- par()$usr
dev.off()
pdf("SiteMapA.pdf", useDingbats=FALSE, width=11, height=10)
  plot(GulfStates_modified, lty=0, col=grey(0.95), xlim=usr[1:2], ylim=usr[3:4])
  plot(spRESTORE_MGCV_BND, col="white", add=TRUE, lty=0)
  plot(edwards_aquifer_outcrop, col="#F99B00", lty=0, add=TRUE)
  plot(StreamsOutRESTORE, add=TRUE, lwd=0.15, col="#91B0BD")
  plot(StreamsInRESTORE, add=TRUE, lwd=0.17, col="#6AC3F2")
  plot(GulfStates_modified, add=TRUE, lwd=.4, lty=2, col=NA)
  ix <- length(2:(length(bnd_poly_aea[,1])-1))
  ix <- c(1,sort(sample(ix, size=20000, replace=FALSE)),length(bnd_poly_aea[,1]))
  lines(bnd_poly_aea[ix,1], bnd_poly_aea[ix,2], lty=1, lwd=1.5, col="#006F41")
  tmp <- DD[DD$ed_rch_zone == 1,]
  for(site in unique(tmp$site_no)) {
     tmp2 <- tmp[tmp$site_no == site, ]; tmp2 <- tmp2[1,]
     plot(tmp2, pch=1, lwd=0.7, cex=0.8, col="#8D4200", add=TRUE)
  }
  for(site in unique(DD$site_no)) {
    tmp2 <- DD[DD$site_no == site, ]; tmp2 <- tmp2[1,]
    plot(tmp2, pch=2, lwd=0.7, cex=0.8, col="#1E4D2B", add=TRUE)
  }
  map_annotation()
dev.off()
pdf("SiteMapB.pdf", useDingbats=FALSE, width=11, height=10)
  plot(GulfStates_modified, lty=0, col=grey(0.95), xlim=usr[1:2], ylim=usr[3:4])
  plot(spRESTORE_MGCV_BND, col="white", add=TRUE, lty=0)
  #plot(edwards_aquifer_outcrop, col="#F99B00", lty=0, add=TRUE)
  plot(GulfStates_modified, add=TRUE, lwd=.4, lty=2, col=NA)
  k <- 0; ks <- c(0.6,0.8,1.0,1.2,1.4,1.6)
  for(d in sort(unique(DD$decade))) {
    k <- k + 1
    plot(DD[DD$decade == d,], pch=1, lwd=0.7, cex=ks[k], col=rgb(1,0,0,.5), add=TRUE)
  }
  ix <- length(2:(length(bnd_poly_aea[,1])-1))
  ix <- c(1,sort(sample(ix, size=20000, replace=FALSE)),length(bnd_poly_aea[,1]))
  lines(bnd_poly_aea[ix,1], bnd_poly_aea[ix,2], lty=1, lwd=1.5, col="#006F41")
  plot(DD, pch=2, lwd=0.7, cex=0.6, col="#1E4D2B", add=TRUE)
  map_annotation()
dev.off()

