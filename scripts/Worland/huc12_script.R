library(jsonlite)
library(dplyr)
setwd("~/Documents/Projects/WaterSmart/5_data/gages_3/")
huc_df<-readr::read_csv(file = "huc12/huc12_comid.csv")
comids<-huc_df[2]
names(comids) <- "COMID"
metadata<-fromJSON("../databaseShapefiles/catchment_characteristics/data_cleanup/metadata.json")
datasets<-list()
extension<-"_tot.rds"
for(url in unique(metadata$datasteURL)) {
  theme<-metadata$themeLabel[min(which(metadata$datasteURL %in% url))]
  rdsFile<-paste0("../databaseShapefiles/catchment_characteristics/data_cleanup/rds/broken_out/",
                  strsplit(url,split = "/")[[1]][6],extension)
  outFile<-paste0("huc12/broken_out/huc12_",strsplit(url,split = "/")[[1]][6],extension)
  if(!"Chemical" %in% theme) {
    if(!file.exists(outFile)) {
      print(paste("Creating subset of", metadata$datasetLabel[min(which(metadata$datasteURL %in% url))]))
      varData<-readRDS(rdsFile)
      if(length(names(varData))>1) {
        inData<-left_join(comids, varData, by = c("COMID"))
        saveRDS(inData, outFile)
      } else {
          print("Didn't find data for this dataset.")
      }
    } else {
        print(metadata$datasetLabel[min(which(metadata$datasteURL %in% url))])
        inData<-readRDS(outFile)
    }
    if(length(names(inData))>1) {
      huc_df<-cbind(huc_df,inData[2:length(names(inData))])
    } else {
      print("Didn't find data for this dataset.")
    }
  }
}
saveRDS(huc_df, file=paste0("huc12/huc12_mw_characteristics",extension))