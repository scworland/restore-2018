library(sp)
library(rgdal)

ALBEA <- paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 ",
     "+x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
ALBEA <- sp::CRS(ALBEA)

unzip("restore_mgcv_bnd.zip", overwrite=TRUE)
huc_outline <- readOGR("restore_mgcv_bnd", "restore_mgcv_bnd")
unlink("./restore_mgcv_bnd/", recursive=TRUE)

huc_outline <- spTransform(huc_outline, ALBEA)

# one role of this script that is no longer needed is the isolation of the
# parent polygon. Original processing showed that tiny errors in the geometry
# exists by micro islands of polygons.
# These were removed and the restore_mgcv_bnd resaved.

h <- huc_outline
hp <- slot(h, "polygons")
hh <- hp[[1]]
jj <- slot(hh, "Polygons")
gg <- lapply(jj, function(i) slot(i, "coords"))
lb <- lapply(jj, function(i) { slot(i, "labpt") })
for(i in 1:length(lb)) { tmp <- lb[[i]]; text(tmp[1], tmp[2], i, col=2) }
bb <- slot(jj[[1]], "area")
xy <- slot(jj[[1]], "coords")
x <- xy[,1]; y <- xy[,2]
soap <- SpatialPolygons(list(Polygons(list(Polygon(xy)), ID = 1)), proj4string = ALBEA)

spRESTORE_MGCV_BND <- SpatialPolygonsDataFrame(soap,
                                               data.frame(type="RESTORE BOUNDARY FOR MGCV(gam)") )

proj4string(spRESTORE_MGCV_BND) <- ALBEA
#writeOGR(spRESTORE_MGCV_BND, "restore_mgcv_bnd/", "restore_mgcv_bnd",
#                             driver="ESRI Shapefile")

xy <- coordinates(spRESTORE_MGCV_BND)
bnd_poly_aea <- xy
save(spRESTORE_MGCV_BND, bnd_poly_aea, file="RESTORE_MGCV_BND.RData")


unzip("restore_bnd.zip", overwrite=TRUE)
spBND <- readOGR("restore_bnd", "restore_bnd")
unlink("./restore_bnd/", recursive=TRUE)
save(spBND, file="RESTORE_BND.RData")
if(dir.exists("./__MACOSX")) unlink("./__MACOSX", recursive=TRUE)
