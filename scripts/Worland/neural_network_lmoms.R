
library(keras)
library(tidyverse)
library(feather)
library(rethinking)

# setwd("~/Documents/Restore")

d <- read_feather("data/gage/all_gage_data.feather")

Y <- select(d,pplo,L1:T4) %>% 
  mutate(T2 = L2/L1,
         L1 = log10(L1)) %>%
  select(L1,T2,T3,T4,pplo) 

X <- select(d,ppt_mean:statsgo) %>%
  mutate_at(vars(bedperm:statsgo),funs(as.numeric(as.factor(.)))) %>%
  mutate_all(funs(as.vector(scale(.)))) %>%
  select_if(~!any(is.na(.))) %>% 
  as.matrix()

# Build model function for multiple outputs
build_model <- function(){
  input <- layer_input(shape=dim(X)[2],name="basinchars")
  
  base_model <- input  %>%
    layer_dense(units = 32, activation='relu') %>%
    layer_dropout(rate=0.0) 

  l1_pred <- base_model %>% 
    layer_dense(units = 1, activation='relu', name="l1") 
  
  t2_pred <- base_model %>% 
    layer_dense(units = 1, activation='relu', name="t2") 
  
  t3_pred <- base_model %>% 
    layer_dense(units = 1, activation='relu', name="t3") 
  
  t4_pred <- base_model %>% 
    layer_dense(units = 1, activation='relu', name="t4") 
  
  pplo_pred <- base_model %>% 
    layer_dense(units = 1, activation='softplus', name="pplo") 
  
  model <- keras_model(input,list(l1_pred,t2_pred,t3_pred,t4_pred,pplo_pred)) %>%
    compile(optimizer = "rmsprop",
            loss="mse",
            metrics="mae",
            loss_weights=c(0.2,0.5,0.5,1,1))
  
  return(model)
}

# K-fold cross validation
source("scripts/Worland/keras_funs.R")

cv_results <- k_foldcv(k=10,epochs=125,batch_size=50,Y=Y,X=X,data=d)

tail(cv_results$average_mse)

ggplot(cv_results$average_mse) +
  geom_line(aes(x = epoch, y = val_mse))

# plots estimated vs observed
est_obs <- cv_results$est_obs %>%
  # mutate(nnet = ifelse(variable=="L1",10^nnet,nnet),
  #        obs = ifelse(variable=="L1",10^obs,obs)) %>%
  gather(model,value,-variable,-obs,-site_no,-decade) %>%
  mutate(col=ifelse(obs>value,1,0))

ggplot(est_obs) +
  geom_point(aes(obs,value,color=col),alpha=0.2) +
  geom_abline(slope=1,intercept=0,linetype="dashed",color="black") +
  facet_wrap(~variable,scales = "free",ncol=2) +
  labs(x="observed L-moment",y="estimated L-moment") +
  ggtitle("10 fold CV") +
  theme_bw() +
  theme(legend.position = 'none')

# subtitle="where drainage area = mean(log10(lmom[-k])/log10(area[-k]))*log10(area[k])"

est_obs <- cv_results$est_obs %>%
  mutate(nnet = ifelse(variable=="L1",10^nnet,nnet),
         obs = ifelse(variable=="L1",10^obs,obs)) %>%
  select(-obs) %>%
  spread(variable,nnet) %>%
  select(site_no,decade,L1,T2,T3,T4,pplo)

write_feather(est_obs,"data/gage/lmom_estimates.feather")
  