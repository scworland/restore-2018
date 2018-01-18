
library(keras)
library(tidyverse)
library(feather)
library(rethinking)
library(forcats)

setwd("~/Documents/Restore")

d <- read_feather("data/gage/all_gage_data.feather")

Y <- select(d,L1:T3) %>% 
  mutate(L1 = log10(L1),
         L2 = log10(L2),
         T3 = T3) %>%
  as.matrix()

X <- select(d,ppt_mean:tot_rdx) %>%
  mutate_all(funs(as.vector(scale(.)))) %>%
  select_if(~!any(is.na(.))) %>% 
  as.matrix()

# Build model function for multiple outputs
build_model <- function(){
  input <- layer_input(shape=dim(X)[2],name="basinchars")
  
  base_model <- input  %>%
    layer_dense(units = 64, activation='relu') %>%
    layer_dropout(rate=0.1) %>%
    layer_dense(units = 64, activation='relu') %>%
    layer_dropout(rate=0.1)
  
  l1_pred <- base_model %>% 
    layer_dense(units = 1, name="l1") 
  
  l2_pred <- base_model %>% 
    layer_dense(units = 1, name="l2") 
  
  t3_pred <- base_model %>% 
    layer_dense(units = 1, name="t3") 
  
  model <- keras_model(input,list(l1_pred,l2_pred,t3_pred)) %>%
    compile(optimizer = "rmsprop",
            loss="mse",
            metrics="mae",
            loss_weights=c(0.8,0.8,1))
  
  return(model)
}

# K-fold cross validation
source("scripts/Worland/keras_funs.R")

run <- list()
for(i in 1:20) {
  print(paste0(i, " of 20 runs"))
  cv_results <- k_foldcv(k=10,N=nrow(d),epochs=100,batch_size=250,data=d)
  run[[i]] <- select(cv_results$est_obs,-drainage_area)
}

saveRDS(run, "data/sandbox/dropout_20.rds")

run <- readRDS("data/sandbox/dropout_20.rds")

runs_combined <- run %>%
  reduce(left_join, by = c("site_no","decade","lmom","obs"))

run_df <- runs_combined %>%
  select(-site_no,-decade,-obs,-lmom) %>%
  setNames(paste0("run_",1:length(run)))

mu <- apply(run_df,1,mean)
low <- apply(run_df,1,PI,prob=0.9)[1,]
high <- apply(run_df,1,PI,prob=0.9)[2,]

lmoms <- runs_combined %>%
  select(site_no,decade,lmom,obs) %>%
  mutate(mu=mu,
         low=low,
         high=high) %>%
  rownames_to_column(var = "id") %>%
  arrange(lmom) 

ggplot(lmoms[9000:9100,]) +
  geom_pointrange(aes(x=fct_reorder(factor(id),obs),y=mu,ymin=low,ymax=high),alpha=1) +
  geom_point(aes(x=id,y=obs),fill="dodgerblue",color="black",shape=23) +
  ggtitle("Predicted posterior L1 vs observed for 100 sites",
          subtitle="A dropout rate of 0.1 was used to estimate the posterior distribution") +
  ylab("predicted mean and 90% PI") +
  xlab("site") +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
