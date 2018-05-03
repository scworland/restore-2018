
sw_add_housedensity <- function(item_list,sites,huc12s){
  
  filter <- dplyr::filter
  select <- dplyr::select
  sites <- select(sites,-huc12,-lon,-lat)
  huc12s <- select(huc12s,-lon,-lat)
  
  # housing density 
  item_list <- item_list %>%
    filter(item=="5910de31e4b0e541a03ac983") 
  
  houses <- sw_sb_extract(item_list$item,sites=sites,huc12s=huc12s)
  
  # create gage housing density data frame
  house_gage_df <- as.data.frame(houses$gages) %>%
    select(-matches("cat|tot|nodata")) %>%
    left_join(sites, by="comid") %>%
    distinct(site_no,.keep_all = T) %>%
    select(comid,site_no,everything()) %>%
    gather(variable,value,-comid,-site_no) %>%
    mutate(decade=parse_number(variable),
           decade=recode(decade,
                         '40'=1940,
                         '50'=1950,
                         '60'=1960,
                         '70'=1970,
                         '80'=1980,
                         '90'=1990,
                         '0'=2000,
                         '10'=2010),
           variable=gsub("\\d", "", variable)) %>%
    spread(variable,value) %>%
    filter(decade > 1940) 
  
  #write_feather(house_gage_df,"data/basinchars/nhd_sb/housedens_gage_df.feather")
  
  # create huc housing density data frame
  house_huc12_df <- as.data.frame(houses$hucs) %>%
    select(-matches("cat|tot|nodata")) %>%
    left_join(huc12s, by="comid") %>%
    distinct(huc12,.keep_all = T) %>%
    select(comid,huc12,everything()) %>%
    gather(variable,value,-comid,-huc12) %>%
    mutate(decade=parse_number(variable),
           decade=recode(decade,
                         '40'=1940,
                         '50'=1950,
                         '60'=1960,
                         '70'=1970,
                         '80'=1980,
                         '90'=1990,
                         '0'=2000,
                         '10'=2010),
           variable=gsub("\\d", "", variable)) %>%
    spread(variable,value) %>%
    filter(decade > 1940) 
  
  #write_feather(house_huc12_df,"data/basinchars/nhd_sb/housedens_huc12_df.feather")
  
  house_all <- list(house_gage_df=house_gage_df,
                    house_huc12_df=house_huc12_df)
  
  return(house_all)
}
