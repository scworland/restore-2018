library(sp)
library(rgdal)

unzip("restore_edwards_tsms.zip", overwrite=TRUE)
edwards_tsms          <- readOGR("restore_edwards_tsms/",          "edwards_tsms")
EdwardsAquiferOutcrop <- readOGR("restore_edwards_tsms/", "EdwardsAquiferOutcrop")
unlink("./restore_edwards_tsms/", recursive=TRUE)
if(dir.exists("./__MACOSX")) unlink("./__MACOSX", recursive=TRUE)

ALBEA <- paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 ",
     "+x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
ALBEA <- sp::CRS(ALBEA)

edwards_tsms          <- spTransform(edwards_tsms,          ALBEA)
EdwardsAquiferOutcrop <- spTransform(EdwardsAquiferOutcrop, ALBEA)

save(EdwardsAquiferOutcrop, file="EdwardsAquiferOutcrop.RData")
