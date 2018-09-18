
sw_add_immutable_vars <- function(item_list,sites,huc12s){
  
  filter <- dplyr::filter
  select <- dplyr::select
  sites <- select(sites,-huc12,-lon,-lat)
  huc12s <- select(huc12s,-lon,-lat)
  
  # create immutable basin characteristic matrix
  item_list <- item_list %>% 
    filter(type=="immutable") 
  
  items <- item_list$item
  group <- item_list$grouping_var
  sb_type <- item_list$sb_type
  
  immutable_chars <- list()
  immutable_gages <- list()
  immutable_hucs <- list()
  for (i in 1:length(items)){
    immutable_chars[[i]] <- sw_sb_extract(items[i],type=sb_type[i],group=group[i],sites=sites,huc12s=huc12s)
    immutable_gages[[i]] <- immutable_chars[[i]]$gages
    immutable_hucs[[i]] <- immutable_chars[[i]]$hucs
    print(paste0("completed ",item_list$description[i]," = ",i, "/",length(items)))
  }
  
  
  immutable_gage_df <- as.data.frame(immutable_gages) %>%
    rename(COMID=comid) %>%
    select(-contains("comid",ignore.case = F)) %>%
    rename(comid=COMID) %>%
    left_join(sites, by="comid") %>%
    distinct(site_no, .keep_all = T) %>%
    select(-matches("cat|tot|nodata|acc_s1|x1")) %>%
    select(comid, site_no, everything())
  
  
  immutable_huc12_df <- as.data.frame(immutable_hucs) %>%
    rename(COMID=comid) %>%
    select(-contains("comid",ignore.case = F)) %>%
    rename(comid=COMID) %>%
    left_join(huc12s, by="comid") %>%
    distinct(huc12, .keep_all = T) %>%
    select(-matches("cat|tot|nodata|acc_s1|x1")) %>%
    select(comid, huc12, everything())
  
  
  immutable_all <- list(immutable_gage_df=immutable_gage_df,
                        immutable_huc12_df=immutable_huc12_df)
  
  # need a change
  hold <- NULL
  
  return(immutable_all)
  
}
