library(tidyverse)
library(feather)
library(corrr)

source("scripts/Worland/utils.R")

# load data
d <- remake::fetch('gage_all')

Y <- select(d,f50) %>%
  mutate_all(funs(.+2)) %>%
  mutate_all(funs(log10))

X <- select(d,lon,lat,acc_hdens:statsgo,-perennial_ice_snow,-area_sqkm) %>%
  mutate_at(vars(bedperm:statsgo),funs(as.numeric(as.factor(.)))) %>%
  mutate_all(funs(as.vector(scale(.)))) %>%
  select_if(~!any(is.na(.)))

# correlation analysis with quantile
corrr_analysis <- X %>%
  mutate(y = Y$f50) %>%
  correlate() %>%
  focus(y) %>%
  rename(feature = rowname) %>%
  arrange(abs(y)) %>%
  mutate(feature = as.factor(feature))

corrr_analysis %>%
  ggplot(aes(x = y, y = fct_reorder(feature, desc(y)))) +
  geom_point() +
  geom_segment(aes(xend = 0, yend = feature), 
               color = palette_light()[[2]], 
               data = corrr_analysis %>% filter(y > 0)) +
  geom_point(color = palette_light()[[2]], 
             data = corrr_analysis %>% filter(y > 0)) +
  geom_segment(aes(xend = 0, yend = feature), 
               color = palette_light()[[1]], 
               data = corrr_analysis %>% filter(y < 0)) +
  geom_point(color = palette_light()[[1]], 
             data = corrr_analysis %>% filter(y < 0)) +
  labs(x="correlation with 50th percentile",y="features") +
  theme_bw()

# correlation EDA
x_cor <- X %>%
  #select(barren:woody_wetland) %>%
  correlate() %>% 
  shave() %>%
  stretch() %>%
  arrange(desc(abs(r))) %>%
  na.omit()

rplot(x_cor)
