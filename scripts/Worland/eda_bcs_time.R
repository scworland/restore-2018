
library(tidyverse)

setwd("~/Documents/Restore")

# extract LULC file names
list_file <- list.files(path="data/basinchars/",pattern = glob2rx("*LULC*.feather"))

# average basins for each year ----

lulc_comb <- NULL
for (i in 1:length(list_file)){
  fname=list_file[i] # grab first file name 
  options(warn=-1) # turn off warnings
  lulc_hold <- read_feather(paste0("data/basinchars/",fname)) %>% 
    mutate_all(funs(replace(., . == -9999, NA))) %>% # replace NAs
    select(-SITE_NO) %>% # remove siteno for plots
    mutate(year=parse_number(fname)) %>% # extract year from fname
    summarize_all(funs(mean),na.rm=T) # take mean across basins
  lulc_comb <- rbind(lulc_comb,lulc_hold) # append dataframes
  options(warn=0) # turn on warnings
}

# make plot
lulc_comb %>%
  gather(key,value,-year) %>%
  ggplot() + 
  geom_line(aes(year,value)) + 
  facet_wrap(~key, scale="free")

# spot check random basins ----

lulc_comb <- NULL
for (i in 1:length(list_file)){
  fname=list_file[i] # grab first file name 
  options(warn=-1) # turn off warnings
  lulc_hold <- read_feather(paste0("data/basinchars/",fname)) %>% 
    mutate_all(funs(replace(., . == -9999, NA))) %>% # replace NAs
    mutate(year=parse_number(fname)) # extract year from fname

  lulc_comb <- rbind(lulc_comb,lulc_hold) # append dataframes
  options(warn=0) # turn on warnings
}

lulc_rand <- lulc_comb %>%
  mutate(siteno=as.factor(SITE_NO)) %>%
  select(-SITE_NO) %>%
  filter(siteno %in% sample(levels(siteno),5)) %>%
  arrange(siteno) %>%
  gather(key,value,-year,-siteno) %>%
  ggplot() + 
  geom_line(aes(year,value)) + 
  facet_wrap(~key, scale="free")




