

add_huc12_comid <- function(restore_hucs){
  
  # get NHD comids for each huc 12
  huc12_list <- data.frame(HUC_12=restore_hucs$HUC_12) %>%
    mutate(HUC_12 = as.character(HUC_12))
  
  # use NLDI to find comids
  huc12_list$comid <- NA
  for(i in 1:nrow(huc12_list)){
    d <- query_nldi("huc12pp",huc12_list$HUC_12[i],tier = "test")
    if(!is.null(d)) {
      huc12_list$comid[i] <- d$features$properties$comid[1]
    } else {
      print(paste0("no data for ", huc12_list$HUC_12[i]))
    }
    print(paste0("huc #: ", i))
  }
  
  huc12_list <- huc12_list %>%
    mutate(comid = ifelse(comid == "",NA,comid)) %>%
    na.omit()
  
  return(huc12_list)
}