
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
    Q <- qlmomco(FF,par)
  }
  return(Q)
}

sw_left_na <- function(x,num){
  n <- num - length(x)
  if(n > 0){
    result <- c(rep(NA,n),x)
  }else{
    result <- x
  }
  return(result)
}

# find column name with max value
find_max <- function(df) {
  colnames(df)[max.col(df,ties.method="first")]
}

# download and unzip files from sciencebase
#' @example sw_sb_extract("5835cad4e4b0d9329c801b22")
sw_sb_extract <- function(item,type="ACC",path="data/basinchars/nhd_sb",group=0,sites,huc12s){
  
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
    
    if(!file.exists(pth)){
    dat <- item_file_download(item, 
                              names=files$fname[i], 
                              destinations=pth,
                              overwrite_file = T)
    }else{
      dat <- pth
    }
    
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
sw_find_index <- function(all_x,index_x,site,distance=200){
  
  library(proxy)
  
  decades <- unique(all_x$decade)
  
  index_site <- NULL
  subs <- NULL
  for(i in 1:length(decades)){
    
    # extract rows for site (simulate ungaged)
    site_data <- filter(all_x,comid==site & decade==decades[i])
    
    # subset possible reference by decade
    index_xd <- filter(index_x,decade==decades[i])
    
    # possible donor sites
    index_sites <- index_xd$comid
    lon <- index_xd$lon
    lat <- index_xd$lat
    
    # find sites with streamflow and within distance
    neighbors <- sw_geo_dist(site_data,index_sites,lon,lat,distance=distance)
    
    # find site that is the most similar in Euclidean space
    sub_xd <- index_xd %>%
      filter(comid %in% neighbors$neighbor) %>%
      filter(comid != site)
    
    dists <- proxy::dist(select(sub_xd,lon,lat,acc_hdens:acc_rdx),
                         select(site_data,lon,lat,acc_hdens:acc_rdx))
    
    index_site[i] <- sub_xd$comid[which(dists==min(dists))]
    
    sub_comid <- sub_xd %>%
      select(possible_sites=comid,decade) %>%
      mutate(index_site=index_site[i],
             site=site)
    
    subs <- rbind(subs,sub_comid)
  }
  
  index_sites <- data.frame(decade=decades,index_site,stringsAsFactors = F)
  possible_index <- select(subs,site,decade,possible_sites,index_site)
  result <- list(index_sites=index_sites,possible_index=possible_index)
}

# find sites within geographic distance
sw_geo_dist <- function(site_data,index_sites,lon,lat,distance=200){
  
  library(geosphere)
  
  site_coords <- c(site_data$lon,site_data$lat)
  
  coord_mat <- cbind(lon,lat)
  
  dmat <- distm(coord_mat,site_coords,fun = distHaversine)/1000
  
  neighbors <- data.frame(site=site_data$comid,
                          neighbor=index_sites[which(dmat <= distance)],
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

# return rows in prediction matrix that are outside the range
# of the training dataset
sw_within_range <- function(train,predict) {
  
  if(any(lapply(train, class) %in% "factor")) {
    stop("Columns of 'train' dataframe must be either numeric or character class")
  }
  
  if(any(lapply(predict, class) %in% "factor")) {
    stop("Columns of 'predict' dataframe must be either numeric or character class")
  }
  
  index_all <- NULL
  for(i in 1:ncol(predict)){
    name <- names(predict)[i]
    p <- predict[,name]
    t <- train[,name]
    
    if(class(t)=="character"){
      rows <- which(!p %in% t)
      if(!purrr::is_empty(rows)){
        index <- data.frame(variable=name,
                            index=rows,
                            value=p[rows],
                            stringsAsFactors = F)
        
        index_all = rbind(index_all,index)
        
      }else{
        index_all = index_all
      }
    }else{
      rows <- which(!p >= range(t)[1] | !p <= range(t)[2])
      if(!purrr::is_empty(rows)){
        index <- data.frame(variable=name,
                            index=rows,
                            value=p[rows],
                            stringsAsFactors = F)
        index_all = rbind(index_all,index)
      }else{
        index_all = index_all
      }
    }
    
    index_all = rbind(index_all,index)
    
  }
  
  return(index_all)
}

# test function
# train <- data.frame(a = c("a","b","c","d","e"),
#                     b = c(1,2,3,4,5),
#                     c = c(10,11,12,13,14),
#                     stringsAsFactors = F)
# 
# predict <- data.frame(a = c("a","m","c","d","z","b","e"),
#                       b = c(1,2,3,4,5,1.2,6),
#                       c = c(10,11,12,100,14,13,11.5),
#                       stringsAsFactors = F)
# 
# sw_within_range(train,predict)









