
library(keras)
library(tidyverse)
library(feather)
library(rethinking)

# setwd("~/Documents/Restore")

d <- read_feather("data/gage/all_gage_data.feather")

# Y <- select(d,L1:T3) %>% 
#   mutate(L1 = log10(L1),
#          L2 = log10(L2),
#          T3 = T3) %>%
#   as.matrix()

Y <- select(d,L1:T4) %>% 
  mutate(T2 = L2/L1,
         L1 = log10(L1)) %>%
  select(L1,T2,T3,T4) %>%
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
    layer_dense(units = 32, activation='relu') %>%
    layer_dropout(rate=0.1) 
  
  l1_pred <- base_model %>% 
    layer_dense(units = 1, name="l1") 
  
  t2_pred <- base_model %>% 
    layer_dense(units = 1, name="t2") 
  
  t3_pred <- base_model %>% 
    layer_dense(units = 1, name="t3") 
  
  t4_pred <- base_model %>% 
    layer_dense(units = 1, name="t4") 
  
  model <- keras_model(input,list(l1_pred,t2_pred,t3_pred,t4_pred)) %>%
    compile(optimizer = "rmsprop",
            loss="mse",
            metrics="mae",
            loss_weights=c(0.6,0.7,0.7,1))
  
  return(model)
}

# K-fold cross validation
source("scripts/Worland/keras_funs.R")

cv_results <- k_foldcv(k=10,N=nrow(d),epochs=100,batch_size=250,data=d)

tail(cv_results$average_mse)

ggplot(cv_results$average_mse) +
  geom_line(aes(x = epoch, y = val_mse))

# plots estimated vs observed
est_obs <- cv_results$est_obs %>%
  gather(model,value,-lmom,-obs,-site_no,-decade) %>%
  mutate(col=ifelse(obs>value,1,0))

ggplot(est_obs) +
  geom_point(aes(obs,value,color=col),alpha=0.2) +
  geom_abline(slope=1,intercept=0,linetype="dashed",color="black") +
  facet_wrap(model~lmom,scales = "free",ncol=4) +
  labs(x="observed L-moment",y="estimated L-moment") +
  ggtitle("10-fold cross validated predictions",
         subtitle="where drainage area = mean(log10(lmom[-k])/log10(area[-k]))*log10(area[k])") +
  theme_bw() +
  theme(legend.position = 'none')


