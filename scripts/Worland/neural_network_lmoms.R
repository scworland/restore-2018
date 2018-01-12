
library(keras)
library(tidyverse)
library(feather)

d <- read_feather("data/gage/all_gage_data.feather")

Y <- select(d,L1:T3) %>% 
  mutate(L1 = log10(L1),
         L2 = log10(L2)) %>%
  as.matrix()

X <- select(d,ppt_mean:tot_rdx) %>%
  mutate_all(funs(as.vector(scale(.)))) %>%
  select_if(~!any(is.na(.))) %>% 
  as.matrix()

# create train and test set
set.seed(3)
ind <- sample(2, nrow(Y), replace=T, prob=c(0.80, 0.20))

Xtrain <- X[ind==1,]
Ytrain <- Y[ind==1,]
Xtest <- X[ind==2,]
Ytest <- Y[ind==2,]

# Build model for multiple outputs
input <- layer_input(shape=dim(Xtrain)[2],name="basinchars")

base_model <- input  %>%
  layer_dense(units = 64, activation='relu') %>% 
  layer_dense(units = 64, activation='relu') %>%
  layer_dense(units = 64, activation='relu') 

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

model_fit <- model %>% 
  fit(x=Xtrain,
      y=list(Ytrain[,1],Ytrain[,2],Ytrain[,3]), 
      epochs=130, 
      batch_size = 5,
      validation_split = 0.1, 
      verbose=0)

plot(model_fit,metrics=c("loss","l1_loss","l2_loss","t3_loss"),theme_bw=F)

# make predictions
yhat <- predict(model, Xtest) %>%
  data.frame() %>%
  setNames(c("L1","L2","T3")) %>%
  gather(lmom,est) %>%
  mutate(obs = gather(data.frame(Ytest),lmom,obs)$obs) 

ggplot(yhat) + 
  geom_point(aes(obs,est),alpha=0.2) + 
  geom_abline(slope=1,intercept=0,linetype="dashed",color="dodgerblue") +
  facet_wrap(~lmom,scales="free") +
  labs(x="observed L-moment",y="estimated L-moment") +
  theme_bw()



