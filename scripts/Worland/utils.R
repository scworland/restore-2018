
# remove aux sites from dv list
remove_aux <- function(l){
  ncols <- data.frame(ncols=unlist(lapply(l,function(x){ncol(x)})))
  ncols$siteno <- rownames(ncols)
  aux_siteno <- ncols$siteno[which(ncols$ncols !=10)]
  
  for(i in 1:length(aux_siteno)){
    hold <- l[[aux_siteno[i]]]
    hold <- hold[ , -which(names(hold) %in% c("AUX.GAGE_Flow","AUX.GAGE_Flow_cd"))]
    l[[aux_siteno[i]]] <- hold
  }
  
  return(l)
}

# calculate empirical fdcs
sw_efdc <- function(Q){
  df <- data.frame(Q) %>%
    mutate(n = nrow(.),
           m = min_rank(Q),
           p = m/(n+1)) 
  
  p <- as.numeric(df$p)
}

# generate ep over range based on L-moments
sw_lmoms <- function(x, type="pe3") { 
  library(lmomco)
  moms <- lmomco::lmoms(x, no.stop=TRUE)
  if(!are.lmom.valid(moms)) return(NULL) 
  pars <- lmomco::lmom2par(moms, type=type)
  if(is.null(pars)) return(NULL) 
  return(1-plmomco(x, pars))
}

# generate quantiles over range based on L-moments
sw_lmoms2 <- function(x,type="gno") { 
  library(lmomco)
  moms <- vec2lmom(as.numeric(x))
  if(!are.lmom.valid(moms)) return(NULL) 
  pars <- lmomco::lmom2par(moms, type=type)
  if(is.null(pars)) return(NULL) 
  FF <- nonexceeds()
  return(qlmomco(FF, pars))
}

sw_lmom2q <- function(lmoms,FF){
  moms <- vec2lmom(lmoms, lscale=FALSE) 
  par <- paraep4(moms, snap.tau4=TRUE) 
  if(is.null(par)){
    Q <- NA
  }else{
    Q <- qlmomco(FF, par)
  }
  return(Q)
}

# find column name with max value
find_max <- function(df) {
  colnames(df)[max.col(df,ties.method="first")]
}

# download and unzip files from sciencebase
#' @example sw_sb_extract("5835cad4e4b0d9329c801b22")
sw_sb_extract <- function(item,type="TOT",path="data/basinchars/nhd_sb",group=0,sites,huc12s){
  
  library(sbtools)
  
  fnames <- item_list_files(item)$fname
  
  if(any(grepl(type,fnames))) {
    files <- item_list_files(item) %>%
      filter(grepl(type,fname))
  } else {
    files <- item_list_files(item) 
  }
  
  files <- filter(files,!grepl("xml",fname))
  
  gages <- select(sites,comid)
  hucs <- select(huc12s,comid)
  for (i in 1:nrow(files)){
    
    pth <- file.path(path,files$fname[i])
    
    dat <- item_file_download(item, 
                              names=files$fname[i], 
                              destinations=pth,
                              overwrite_file = T)
    
    if(group==1){
      
      d <- read_delim(unzip(dat,exdir=path,overwrite=T),",",guess_max = 20000) %>%
        rename_all(tolower) %>%
        mutate(comid = as.character(comid)) %>%
        mutate_at(vars(-comid), funs(as.numeric)) %>%
        mutate(dominant_group=find_max(.[,-1])) %>%
        select(comid,dominant_group)
      
      colnames(d)[2] <- str_extract(files$fname[1], ".+?(?=_)")
      
    }else{
    
    d <- read_delim(unzip(dat,exdir=path,overwrite=T),",",guess_max = 20000) %>%
      rename_all(tolower) %>%
      mutate(comid = as.character(comid)) %>%
      mutate_at(vars(-comid), funs(as.numeric))
    
    }
    
    gages <- left_join(gages,d,by="comid")
    hucs <-  left_join(hucs,d,by="comid")

  }
  
  out <- list(gages=gages,hucs=hucs)
  return(out)
}

# return sites in same kcluster as given site
sw_obs_in_cluster <- function(site,kmeans_df) {
  
  site_df <- kmeans_df %>%
    filter(site_no==site)
  
  n <- nrow(site_df)
  
  datalist <- list()
  for (i in 1:n){
    df  <- kmeans_df %>%
      select(site_no,decade,kclust) %>%
      filter(kclust==site_df$kclust[i] & decade==site_df$decade[i])
    
    datalist[[i]] <- df
  }
  
  out_df <- data.frame(dplyr::bind_rows(datalist)) %>%
    select(site_no,decade)
}

# Find reference site
sw_find_ref <- function(all_x,ref_x,site,distance=200){
  
  library(proxy)
  
  decades <- unique(all_x$decade)
  
  ref_site <- NULL
  subs <- NULL
  for(i in 1:length(decades)){
    
    # extract rows for site (simulate ungaged)
    site_data <- filter(all_x,site_no==site & decade==decades[i])
    
    # subset possible reference by decade
    ref_xd <- filter(ref_x,decade==decades[i])
    
    # possible donor sites
    ref_sites <- ref_xd$site_no
    lon <- ref_xd$lon
    lat <- ref_xd$lat
    
    # find sites with streamflow and within distance
    neighbors <- sw_geo_dist(site_data,ref_sites,lon,lat,distance=distance)
    
    # find site that is the most similar in Euclidean space
    sub_xd <- ref_xd %>%
      filter(site_no %in% neighbors$neighbor) %>%
      filter(site_no != site)
    
    dists <- proxy::dist(select(sub_xd,ppt_mean:tot_rdx),
                         select(site_data,ppt_mean:tot_rdx))
    
    ref_site[i] <- sub_xd$site_no[which(dists==min(dists))]
    
    sub_site_no <- sub_xd %>%
      select(possible_sites=site_no,decade) %>%
      mutate(ref_site=ref_site[i],
             site=site)
    
    subs <- rbind(subs,sub_site_no)
  }
  
  refs <- data.frame(decade=decades,ref_site,stringsAsFactors = F)
  subs <- select(subs,site,decade,possible_sites,ref_site)
  result <- list(refs=refs,subs=subs)
}

# find sites within geographic distance
sw_geo_dist <- function(site_data,ref_sites,lon,lat,distance=200){
  
  library(geosphere)
  
  site_coords <- c(site_data$lon,site_data$lat)
  
  coord_mat <- cbind(lon,lat)
  
  dmat <- distm(coord_mat,site_coords,fun = distHaversine)/1000
  
  neighbors <- data.frame(site=site_data$site_no,
                          neighbor=ref_sites[which(dmat <= distance)],
                          stringsAsFactors = F)

  
  return(neighbors)
}


# sample groups. Stole directly from kendonB (https://github.com/tidyverse/dplyr/issues/361)
sample_n_groups = function(tbl, size, replace = FALSE, weight = NULL) {
  grps = tbl %>% 
    groups %>% 
    lapply(as.character) %>% 
    unlist
  
  keep = tbl %>% 
    summarise() %>% 
    ungroup() %>% 
    sample_n(size, replace, weight)
  
  tbl %>% 
    right_join(keep, by=grps) %>% 
    group_by_(.dots = grps)
}










