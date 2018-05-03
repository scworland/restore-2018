
sw_add_lulc <- function(historic_class_link,item_list,sites,huc12s){
  
  filter <- dplyr::filter
  select <- dplyr::select
  sites <- select(sites,-huc12,-lon,-lat)
  huc12s <- select(huc12s,-lon,-lat)
  
  # historic lulc, 1940-1990
  historic_class_link <- historic_class_link %>%
    mutate(variable=sub('.*_', '', variable))
  
  item_list <- item_list %>%
    filter(item=="58cbeef2e4b0849ce97dcd61")
  
  historic_lulc <- sw_sb_extract(item_list$item,sites=sites,huc12s=huc12s)
  
  historic_lulc_gages_df <- as.data.frame(historic_lulc$gages) %>%
    select(-matches("nodata")) %>%
    left_join(sites, by="comid") %>%
    distinct(site_no,.keep_all = T) %>%
    select(comid,site_no,everything()) %>%
    gather(variable,value,-comid,-site_no) %>%
    mutate(decade = sub('.*acc_',"",variable) %>%
             sub('_.*',"",.) %>%
             parse_number(.) %>%
             paste0("19",.),
           variable=sub('.*_', '', variable)) %>%
    left_join(historic_class_link, by="variable") %>%
    filter(lulc_class !="Intentionally left blank",
           lulc_class !="Mining") %>%
    select(-variable) %>%
    distinct(site_no,lulc_class,decade,.keep_all = T) %>%
    select(comid,site_no,decade,variable=lulc_class,value) %>%
    spread(variable,value)
  
  historic_lulc_huc12_df <- as.data.frame(historic_lulc$hucs) %>%
    select(-matches("nodata")) %>%
    left_join(huc12s, by="comid") %>%
    distinct(huc12,.keep_all = T) %>%
    select(comid,huc12,everything()) %>%
    gather(variable,value,-comid,-huc12) %>%
    mutate(decade = sub('.*acc_',"",variable) %>%
             sub('_.*',"",.) %>%
             parse_number(.) %>%
             paste0("19",.),
           variable=sub('.*_', '', variable)) %>%
    left_join(historic_class_link, by="variable") %>%
    filter(lulc_class !="Intentionally left blank",
           lulc_class !="Mining") %>%
    select(-variable) %>%
    distinct(huc12,lulc_class,decade,.keep_all = T) %>%
    select(comid,huc12,decade,variable=lulc_class,value) %>%
    spread(variable,value)
  
  # recent modelled lulc
  recent_lulc <- sw_sb_extract("5a5406bee4b01e7be2308855?groupId=JQUERY-FILE-UPLOAD-ac3590ad-84f6-47f6-86cb-0ad6779fbdbb",sites=sites,huc12s=huc12s)
  
  recent_lulc_gages_df <- as.data.frame(recent_lulc$gages) %>%
    select(-matches("nodata")) %>%
    left_join(sites, by="comid") %>%
    distinct(site_no,.keep_all = T) %>%
    select(comid,site_no,everything()) %>%
    gather(variable,value,-comid,-site_no) %>%
    mutate(decade = sub('.*acc_',"",variable) %>%
             sub('_.*',"",.) %>%
             parse_number(.),
           decade = ifelse(nchar(decade)==2, "1990","2000"),
           variable=sub('.*_', '', variable)) %>%
    left_join(historic_class_link, by="variable") %>%
    filter(lulc_class !="Intentionally left blank",
           lulc_class !="Mining") %>%
    select(-variable) %>%
    distinct(site_no,lulc_class,decade,.keep_all = T) %>%
    select(comid,site_no,decade,variable=lulc_class,value) %>%
    group_by(comid,site_no,decade,variable) %>%
    summarize(value=mean(value)) %>%
    spread(variable,value)
  
  recent_lulc_huc12_df <- as.data.frame(recent_lulc$hucs) %>%
    select(-matches("nodata")) %>%
    left_join(huc12s, by="comid") %>%
    distinct(huc12,.keep_all = T) %>%
    select(comid,huc12,everything()) %>%
    gather(variable,value,-comid,-huc12) %>%
    mutate(decade = sub('.*acc_',"",variable) %>%
             sub('_.*',"",.) %>%
             parse_number(.),
           decade = ifelse(nchar(decade)==2, "1990","2000"),
           variable=sub('.*_', '', variable)) %>%
    left_join(historic_class_link, by="variable") %>%
    filter(lulc_class !="Intentionally left blank",
           lulc_class !="Mining") %>%
    select(-variable) %>%
    distinct(huc12,lulc_class,decade,.keep_all = T) %>%
    select(comid,huc12,decade,variable=lulc_class,value) %>%
    group_by(comid,huc12,decade,variable) %>%
    summarize(value=mean(value)) %>%
    spread(variable,value)
  
  # combine everything
  lulc_gage_df <- bind_rows(historic_lulc_gages_df,recent_lulc_gages_df) %>%
    gather(variable,value,-comid,-site_no,-decade) %>%
    group_by(comid,site_no,decade,variable) %>%
    summarize(value=mean(value)) %>%
    spread(variable,value) %>%
    filter(decade != 1940) %>%
    ungroup() %>%
    mutate(decade=as.numeric(decade))
  
  
  lulc_huc12_df <- bind_rows(historic_lulc_huc12_df,recent_lulc_huc12_df) %>%
    gather(variable,value,-comid,-huc12,-decade) %>%
    group_by(comid,huc12,decade,variable) %>%
    summarize(value=mean(value)) %>%
    spread(variable,value) %>%
    filter(decade != 1940) %>%
    ungroup() %>%
    mutate(decade=as.numeric(decade))
  
  lulc_all <- list(lulc_gage_df=lulc_gage_df,
                   lulc_huc12_df=lulc_huc12_df)
  
  return(lulc_all)
  
}

