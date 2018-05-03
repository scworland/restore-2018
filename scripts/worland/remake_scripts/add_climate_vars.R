

sw_add_climate <- function(item_list,sites,huc12s){
  
  filter <- dplyr::filter
  select <- dplyr::select
  sites <- select(sites,-huc12,-lon,-lat)
  huc12s <- select(huc12s,-lon,-lat)
  
  # temperature ----
  temp_item_list <- item_list %>%
    filter(item=="5787ea72e4b0d27deb377b6d")
  
  temp <- sw_sb_extract(temp_item_list$item,sites=sites,huc12s=huc12s)
  
  # create gage temp data frame
  temp_gage_df <- as.data.frame(temp$gages) %>%
    select(-matches("cat|tot")) %>%
    left_join(sites, by="comid") %>%
    distinct(site_no,.keep_all = T) %>%
    gather(variable,value,-comid,-site_no) %>%
    mutate(year=parse_number(variable),
           decade=as.numeric(sub("\\d$", "", year))*10) %>%
    group_by(comid,site_no,decade) %>%
    summarize(temp_mean = mean(value),
              temp_sd = sd(value)) %>%
    ungroup() %>%
    filter(decade != 1940) %>%
    select(comid, site_no, everything())
  
  # create huc12 temp data frame
  temp_huc12_df <- as.data.frame(temp$hucs) %>%
    select(-matches("cat|tot")) %>%
    left_join(huc12s, by="comid") %>%
    distinct(huc12,.keep_all = T) %>%
    gather(variable,value,-comid,-huc12) %>%
    mutate(year=parse_number(variable),
           decade=as.numeric(sub("\\d$", "", year))*10) %>%
    group_by(comid,huc12,decade) %>%
    summarize(temp_mean = mean(value),
              temp_sd = sd(value)) %>%
    ungroup() %>%
    filter(decade != 1940) %>%
    select(comid, huc12, everything())
  
  # precip ----
  precip_item_list <- item_list %>%
    filter(item=="57bf5c07e4b0f2f0ceb75b1b") 
  
  precip <- sw_sb_extract(precip_item_list$item,sites=sites,huc12s=huc12s)
  
  # create gage precip data frame
  precip_gage_df <- as.data.frame(precip$gages) %>%
    select(-matches("cat|tot")) %>%
    left_join(sites, by="comid") %>%
    distinct(site_no,.keep_all = T) %>%
    gather(variable,value,-comid,-site_no) %>%
    mutate(year=parse_number(variable),
           decade=as.numeric(sub("\\d$", "", year))*10) %>%
    group_by(comid,site_no,decade) %>%
    summarize(ppt_mean = mean(value),
              ppt_sd = sd(value)) %>%
    ungroup() %>%
    filter(decade != 1940) %>%
    select(comid, site_no, everything())

  # create huc12 precip data frame
  precip_huc12_df <- as.data.frame(precip$hucs) %>%
    select(-matches("cat|tot")) %>%
    left_join(huc12s, by="comid") %>%
    distinct(huc12,.keep_all = T) %>%
    gather(variable,value,-comid,-huc12) %>%
    mutate(year=parse_number(variable),
           decade=as.numeric(sub("\\d$", "", year))*10) %>%
    group_by(comid,huc12,decade) %>%
    summarize(ppt_mean = mean(value),
              ppt_sd = sd(value)) %>%
    ungroup() %>%
    filter(decade != 1940) %>%
    select(comid, huc12, everything())
  
  # runoff ----
  runoff_item_list <- item_list %>%
    filter(item=="57bf5e25e4b0f2f0ceb75b77") 
  
  runoff <- sw_sb_extract(runoff_item_list$item,sites=sites,huc12s=huc12s)
  
  # create gage runoff data frame
  runoff_gage_df <- as.data.frame(runoff$gages) %>%
    select(-matches("cat|tot")) %>%
    left_join(sites, by="comid") %>%
    distinct(site_no,.keep_all = T) %>%
    gather(variable,value,-comid,-site_no) %>%
    mutate(year=parse_number(variable),
           decade=as.numeric(sub("\\d$", "", year))*10) %>%
    group_by(comid,site_no,decade) %>%
    summarize(runoff_mean = mean(value),
              runoff_sd = sd(value)) %>%
    ungroup() %>%
    filter(decade != 1940) %>%
    select(comid, site_no, everything())
  
  # create huc12 runoff data frame
  runoff_huc12_df <- as.data.frame(runoff$hucs) %>%
    select(-matches("cat|tot")) %>%
    left_join(huc12s, by="comid") %>%
    distinct(huc12,.keep_all = T) %>%
    gather(variable,value,-comid,-huc12) %>%
    mutate(year=parse_number(variable),
           decade=as.numeric(sub("\\d$", "", year))*10) %>%
    group_by(comid,huc12,decade) %>%
    summarize(runoff_mean = mean(value),
              runoff_sd = sd(value)) %>%
    ungroup() %>%
    filter(decade != 1940) %>%
    select(comid, huc12, everything())
  

  # combine temp, precip, and runoff
  climate_gage_df <- precip_gage_df %>%
    left_join(temp_gage_df,by=c("site_no","comid","decade")) %>%
    left_join(runoff_gage_df,by=c("site_no","comid","decade")) %>%
    distinct(site_no,decade,.keep_all = T) %>%
    select(comid, site_no, everything())

  climate_huc12_df <- precip_huc12_df %>%
    left_join(temp_huc12_df,by=c("huc12","comid","decade")) %>%
    left_join(runoff_huc12_df,by=c("huc12","comid","decade")) %>%
    distinct(huc12,decade,.keep_all = T) %>%
    select(comid, huc12, everything())

  climate_all <- list(climate_gage_df=climate_gage_df,
                      climate_huc12_df=climate_huc12_df)

  return(climate_all)

}