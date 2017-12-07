
library(tidyverse)
library(feather)
library(MultivariateRandomForest)

d <- read_feather("data/gage/lmom3.feather")

Y <- select(d,lm1,lm2,lm3) %>%
  as.matrix()

X <- select(d,siteno,cat_area_sqkm:tot_wildfire_2011) %>%
  gather(variable,value,-siteno) %>%
  group_by(siteno) %>%
  filter(value > 0) %>% 
  ungroup() %>%
  spread(variable, value) %>%
  select_if(~!any(is.na(.))) %>%
  select(-siteno) %>%
  mutate_all(funs(scale)) %>%
  as.matrix()
  
# hold <- CrossValidation(X, Y, 5)

Yhat <- build_forest_predict(X, Y, n_tree=100, m_feature=20, min_leaf=5, X)

output <- data.frame(Y,Yhat)
pars <- lmomco::lmom2par(moms, type="pe3")

