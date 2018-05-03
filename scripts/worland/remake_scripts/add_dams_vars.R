
sw_add_dams <- function(item_list,sites,huc12s){
  
  filter <- dplyr::filter
  select <- dplyr::select
  sites <- select(sites,-huc12,-lon,-lat)
  huc12s <- select(huc12s,-lon,-lat)
  
  # national inventory of damns ----
  item_list <- item_list %>%
    filter(item=="58c301f2e4b0f37a93ed915a") 
  
  dams <- sw_sb_extract(item_list$item,sites=sites,huc12s=huc12s)
  
  # create gage dam data frame using site_no
  dam_gage_df <- as.data.frame(dams$gages) %>%
    select(-matches("cat|tot|dens|.y|.x|2013")) %>%
    left_join(sites, by="comid") %>%
    distinct(site_no,.keep_all = T) %>%
    select(comid,site_no,everything()) %>%
    gather(variable,value,-comid,-site_no) %>%
    mutate(decade=parse_number(variable),
           variable=gsub("\\d", "", variable)) %>%
    spread(variable,value) %>%
    filter(decade > 1940) 
  
  # create create huc12 dam data frame
  dam_huc12_df <- as.data.frame(dams$hucs) %>%
    select(-matches("cat|tot|dens.y|.x|2013")) %>%
    left_join(huc12s, by="comid") %>%
    distinct(huc12,.keep_all = T) %>%
    select(comid,huc12,everything()) %>%
    gather(variable,value,-comid,-huc12) %>%
    mutate(decade=parse_number(variable),
           variable=gsub("\\d", "", variable)) %>%
    spread(variable,value) %>%
    filter(decade > 1940) 
  
  dams_all <- list(dam_gage_df=dam_gage_df, 
                   dam_huc12_df=dam_huc12_df)
  
  return(dams_all)
}
