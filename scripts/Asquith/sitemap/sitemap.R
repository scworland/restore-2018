library(sp)
library(rgeos)
library(rgdal)
library(GISTools)
library(feather)

load("../../../../gisbig/RESTOREstreams.RData") # RESTOREstreams.RData
load("../../../gis/GulfStates.RData") # GulfStates.RData
load("../../../gis/RESTORE_MGCV_BND.RData") # spRESTORE_MGCV_BND.RData
load("../../../gis/EdwardsAquiferOutcrop.RData") # Edwards
load("../fdc_model/FDCEST.RData") # to get bnd
MM <- read_feather("../../../data/gage/all_gage_data.feather") # all_gage_data.feather
MM$ed_rch_zone[MM$site_no == "08155300"] <- 0 # WHA 09/28/2018
MM$ed_rch_zone[MM$site_no == "08155400"] <- 0 # WHA 09/28/2018
MM$ed_rch_zone[MM$site_no == "08156800"] <- 0 # WHA 09/28/2018
MM$ed_rch_zone[MM$site_no == "08181400"] <- 0 # WHA 09/28/2018

LATLONG <- paste0("+init=epsg:4269 +proj=longlat +ellps=GRS80 ",
                  "+datum=NAD83 +no_defs +towgs84=0,0,0")
LATLONG <- sp::CRS(LATLONG)
ALBEA <- paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 ",
                "+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
ALBEA <- sp::CRS(ALBEA)


source("../pred_huc12/a_basemap_funcs.R")


map_base_site_map <- function(xlim=NA, ylim=NA) {
  par(lend=1, ljoin=1)
  polygon(bnd[[1]]$x*1000,bnd[[1]]$y*1000, col=grey(1), lwd=0)
}
  #ix <- length(2:(length(bnd_poly_aea[,1])-1))
  #ix <- c(1,sort(sample(ix, size=20000, replace=FALSE)),length(bnd_poly_aea[,1]))
  #lines(bnd_poly_aea[ix,1], bnd_poly_aea[ix,2], lty=1, lwd=1.5, col="#006F41")


MM <- SpatialPointsDataFrame(cbind(MM$dec_long_va, MM$dec_lat_va), data=MM, proj4string=LATLONG)
MM <- spTransform(MM, ALBEA)


pdf("SiteMap_junk.pdf", useDingbats=FALSE, width=11, height=10)
  plot(spRESTORE_MGCV_BND)
  usr <- par()$usr
dev.off()
unlink("SiteMap_junk.pdf")  # just quietly throw the file away

pdf("SiteMapA.pdf", useDingbats=FALSE, width=11, height=10)
  plot(GulfStates_modified, lty=0, col=grey(0.95), xlim=usr[1:2], ylim=usr[3:4])
  polygon(bnd[[1]]$x*1000,bnd[[1]]$y*1000, col=grey(1), lwd=1.7, lty=0)
  plot(EdwardsAquiferOutcrop, col="#F99B00", lty=0, add=TRUE)
  plot(StreamsOutRESTORE, add=TRUE, lwd=.20, col="#91B0BD")
  plot(StreamsInRESTORE,  add=TRUE, lwd=.25, col="#6AC3F2")
  lines(bnd[[1]]$x*1000,bnd[[1]]$y*1000, col=1, lwd=1.7, lty=1)
  #plot(GulfStates_modified, add=TRUE, lwd=.4, lty=2, col=NA)
  tmp <- MM[MM$ed_rch_zone == 1,]
  for(site in unique(tmp$site_no)) {
     tmp2 <- tmp[tmp$site_no == site, ]; tmp2 <- tmp2[1,]
     plot(tmp2, pch=2, lwd=0.9, cex=0.8, col="#8D4200", add=TRUE)
  }
  for(site in unique(MM$site_no[MM$ed_rch_zone != 1])) {
    tmp2 <- MM[MM$site_no == site, ]; tmp2 <- tmp2[1,]
    plot(tmp2, pch=2, lwd=0.9, cex=0.8, col="#1E4D2B", add=TRUE)
  }
  map_annotation(); legend_est(sitemap=TRUE, triangle=TRUE)
dev.off()
pdf("SiteMapB.pdf", useDingbats=FALSE, width=11, height=10)
  plot(GulfStates_modified, lty=0, col=grey(0.95), xlim=usr[1:2], ylim=usr[3:4])
  polygon(bnd[[1]]$x*1000,bnd[[1]]$y*1000, col=grey(1), lwd=1.7, lty=0)
  plot(EdwardsAquiferOutcrop, col="#F99B00", lty=0, add=TRUE)
  plot(StreamsOutRESTORE, add=TRUE, lwd=.20, col="#91B0BD")
  plot(StreamsInRESTORE,  add=TRUE, lwd=.25, col="#6AC3F2")
  lines(bnd[[1]]$x*1000,bnd[[1]]$y*1000, col=1, lwd=1.7, lty=1)

  #plot(GulfStates_modified, lty=0, col=grey(0.95), xlim=usr[1:2], ylim=usr[3:4])
  #map_base_site_map();
  #plot(GulfStates_modified, add=TRUE, lwd=.4, lty=2, col=NA)
  k <- 0; ks <- c(0.6,0.8,1.0,1.2,1.4,1.6)
  #for(d in sort(unique(DD$decade))) {
  #  k <- k + 1
  #  plot(DD[DD$decade == d,], pch=1, lwd=0.7, cex=ks[k], col="#006F41", add=TRUE)
  #}
  tmp <- D[D$ed_rch_zone != 1,]
  for(d in sort(unique(tmp$decade))) { k <- k + 1
    #tmp2 <- tmp[tmp$site_no == site, ]; tmp2 <- tmp2[1,]
    plot(tmp[tmp$decade == d,], pch=1, lwd=0.7, cex=ks[k], col="#006F41", add=TRUE)
  }
  k <- 0; ks <- c(0.6,0.8,1.0,1.2,1.4,1.6)
  tmp <- DD[DD$ed_rch_zone == 1,]
  for(d in sort(unique(tmp$decade))) { k <- k + 1
     #tmp2 <- tmp[tmp$site_no == site, ]; tmp2 <- tmp2[1,]
     #plot(tmp[tmp$decade == d,], pch=1, lwd=0.7, cex=ks[k], col="#8D4200", add=TRUE)
  }
  map_annotation(); legend_est(sitemap=TRUE, triangle=FALSE, noedwards=TRUE)
  xy <- coordinates(D[D$site_no == "08167000",1])
  points(xy[,1], xy[,2], pch=0, col="#b51b96", cex=2)
  arrows(x0=147500, y0=573000, xy[,1], xy[,2],
         lwd=.5, length=0.15, angle=15, col="#b51b96")
  text(147500, 573000, "Streamgage 08167000 used as an example in the text",
                       cex=0.7, pos=4, col="#7a1265")
dev.off()

cropthem <- TRUE
crop_em <- function(spawn=FALSE) {
  files <- c("SiteMapA", "SiteMapB")
  for(file in files) { my.file <- paste0(file,".pdf")
    system(paste0("pdfcrop --margins '-46 -110 -43 0' --clip ",my.file," ",my.file))
  }
}
crop_em(spawn=cropthem)

print(length(MM$site_no))
print(length(unique(MM$site_no)))


print(length(D$site_no))
print(length(unique(D$site_no)))


