
library(keras)
library(tidyverse)
library(feather)
library(broom)
library(zoo)
source("scripts/Worland/utils.R")

# setwd("~/Documents/Restore")
d <- read_feather("data/gage/all_gage_data2.feather") 

Y <- select(d,f0.02:f99.98) %>%
  mutate_all(funs(.+2)) %>%
  mutate_all(funs(log10))

X <- select(d,lon,lat,acc_hdens:statsgo,-perennial_ice_snow,-area_sqkm) %>%
  mutate_at(vars(bedperm:statsgo),funs(as.numeric(as.factor(.)))) %>%
  mutate_all(funs(as.vector(scale(.)))) %>%
  select_if(~!any(is.na(.))) %>% 
  as.matrix()

# train NN model 
input <- layer_input(shape=dim(X)[2],name="basinchars")

base_model <- input  %>%
  layer_dense(units = 40,activation="relu") %>%
  layer_dropout(rate=0.1) %>%
  layer_dense(units = 30,activation="relu") %>%
  layer_dropout(rate=0.1) 

for(i in 1:dim(Y)[2]){
  y <- colnames(Y)[i]
  outstring <- paste0(
    sprintf("%s <- base_model %%>%%", y), 
    sprintf(" layer_dense(units = 1, activation='relu', name='%s')",y)
  )
  eval(parse(text=outstring))
}

Ylist <- paste0("list(",paste(colnames(Y),sep="",collapse=","),")")
model <- keras_model(input,eval(parse(text=Ylist))) %>%
  compile(optimizer = "rmsprop",
          loss="mse",
          metrics="mae")

model_fit <- model %>% 
  fit(x=X,
      y=Y, 
      epochs=200, 
      batch_size = 30,
      validation_split=0.1, 
      verbose=0)

# data.frame(val_mse=model_fit$metrics$val_loss,epoch=1:150) %>%
#   ggplot() + geom_line(aes(x = epoch, y = val_mse))

# make predictions for every site and decade
covariates <- read_feather("data/gage/all_gage_covariates2.feather") 

all_x <- covariates %>%
  select(lon,lat,acc_hdens:statsgo,-perennial_ice_snow,-area_sqkm,-decade) %>%
  mutate_at(vars(bedperm:statsgo),funs(as.numeric(as.factor(.)))) %>%
  mutate_all(funs(as.vector(scale(.)))) %>%
  select_if(~!any(is.na(.))) %>% 
  as.matrix()

yhat <- predict(model, all_x) %>%
  data.frame() %>%
  setNames(colnames(Y)) %>%
  mutate(site_no=covariates$site_no,
         decade=covariates$decade) %>%
  gather(f,q,-site_no,-decade) %>%
  mutate(q = round((10^q)-2,2),
         f = as.numeric(substring(f, 2))) %>%
  group_by(site_no,decade) %>%
  mutate(q = ifelse(q<0,0,q) %>%
           ifelse(.==cummax(.),.,NA) %>%
           na.approx(.,rule=2)) 

write_feather(yhat,"data/gage/all_fdc_direct_est.feather")
