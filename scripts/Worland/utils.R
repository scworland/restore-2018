
read_months <- function(list,type){
  out <- list()
  for(i in 1:length(list)){
    fname=list[i] # grab first file name 
    month=tolower(substr(fname, start=nchar(fname)-10,stop=nchar(fname)-8))
    options(warn=-1) # turn off warnings
    
    hold <- read_feather(fname) %>%
      select(siteno,contains("CAT")) %>%
      setNames(c("siteno",paste0(type,"_",month)))
    
    out[[i]] <- hold
    
    options(warn=0) # turn on warnings
  }
  return(out)
}

read_years <- function(list){
  out <- list()
  for (i in 1:length(list)){
    fname=list[i] # grab first file name 
    year=parse_number(fname)
    
    options(warn=-1) # turn off warnings
    hold <- read_feather(fname) %>% 
      mutate_all(funs(replace(., . == -9999, NA))) %>% # replace NAs
      rename_at(vars(-contains("SITE")),funs(paste0(.,"_", year))) %>%
      rename(siteno=SITE_NO)
    
    out[[i]] <- hold
    options(warn=0) # turn on warnings
  }
  return(out)
}
