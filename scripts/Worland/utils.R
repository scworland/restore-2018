
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
           m = min_rank(desc(Q)),
           p = m/(n+1)) 
  
  p <- as.numeric(df$p)
}

# generate quantiles over range based on L-moments
sw_lmoms <- function(x, type="pe3") { 
  moms <- lmomco::lmoms(x, no.stop=TRUE)
  if(!are.lmom.valid(moms)) return(NULL) 
  pars <- lmomco::lmom2par(moms, type=type)
  if(is.null(pars)) return(NULL) 
  return(1-plmomco(x, pars))
}

# find column name with max value
find_max <- function(df) {
  colnames(df)[max.col(df,ties.method="first")]
}

# download and unzip files from sciencebase
#' @example sw_sb_extract("5835cad4e4b0d9329c801b22")
sw_sb_extract <- function(item,type="TOT",path="data/basinchars/nhd_sb",group=0){
  
  library(sbtools)
  
  fnames <- item_list_files(item)$fname
  
  if(any(grepl(type,fnames))) {
    files <- item_list_files(item) %>%
      filter(grepl(type,fname))
  } else {
    files <- item_list_files(item) 
  }
  
  files <- filter(files,!grepl("xml",fname))
  
  gages <- select(sites,COMID=comid)
  hucs <- select(huc12s,COMID=comid)
  for (i in 1:nrow(files)){
    
    pth <- file.path(path,files$fname[i])
    
    dat <- item_file_download(item, 
                              names=files$fname[i], 
                              destinations=pth,
                              overwrite_file = T)
    
    if(group==1){
      
      d <- read_delim(unzip(dat,exdir=path,overwrite=T),",",guess_max = 20000) %>%
        mutate(COMID = as.character(COMID)) %>%
        mutate_at(vars(-COMID), funs(as.numeric)) %>%
        mutate(dominant_group=find_max(.[,-1])) %>%
        select(COMID,dominant_group)
      
      colnames(d)[2] <- str_extract(files$fname[1], ".+?(?=_)")
      
    }else{
    
    d <- read_delim(unzip(dat,exdir=path,overwrite=T),",",guess_max = 20000) %>%
      mutate(COMID = as.character(COMID)) %>%
      mutate_at(vars(-COMID), funs(as.numeric))
    
    }
    
    gages <- left_join(gages,d,by="COMID")
    hucs <-  left_join(hucs,d,by="COMID")

  }
  
  out <- list(gages=gages,hucs=hucs)
  return(out)
}

