
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

X <- select(d,lon,lat,acc_hdens:statsgo) %>%
  mutate_at(vars(bedperm:statsgo),funs(as.numeric(as.factor(.)))) %>%
  mutate_all(funs(as.vector(scale(.)))) %>%
  select_if(~!any(is.na(.))) %>% 
  as.matrix()

# Build model function for multiple outputs
build_model <- function(){
  input <- layer_input(shape=dim(X)[2],name="basinchars")
  
  base_model <- input  %>%
    layer_dense(units = 40,activation="relu") %>%
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 30,activation="relu") %>%
    layer_dropout(rate = 0.1)
  
  output <- base_model %>%
    layer_dense(units = 1)
    
  model <- keras_model(input,output) %>%
    compile(optimizer = "rmsprop",
            loss="mse",
            metrics="mae")
  
  return(model)
}

# K-fold cross validation
source("scripts/Worland/keras_funs.R")

preds <- NULL
for(j in 1:ncol(Y)){
  print(paste0("Processing quantile ",j," out of 27"))
  y <- Y[,j]
  cv_results <- k_foldcv(k=10,epochs=200,batch_size=30,Y=y,X=X,data=d)
  preds[[j]] <- cv_results$est_obs
}

preds_df <- bind_rows(preds) %>%
  mutate(obs = (10^obs)-2,
         nnet = round((10^nnet)-2,2),
         variable = as.numeric(substring(variable, 2)))

write_feather(preds_df,"data/gage/direct_indv_estimates.feather")

est_obs <- bind_rows(preds) %>%
  mutate(obs = (10^obs)-2,
         nnet = round((10^nnet)-2,2),
         variable = as.numeric(substring(variable, 2))) %>%
  gather(model,value,-site_no,-decade,-variable) %>%
  group_by(site_no) %>%
  sample_n_groups(6)

ggplot(est_obs) +
  geom_line(aes(variable,value,linetype=model)) +
  scale_linetype_manual(values=c("dashed","solid","dotted")) +
  facet_wrap(site_no~decade,scales="free_y") +
  labs(x="EP",y="Q") +
  theme_bw() +
  scale_y_log10()
