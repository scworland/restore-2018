
library(keras)
library(tidyverse)
library(feather)

d <- read_feather("data/gage/lmom3.feather")

y <- as.matrix(select(d,lm1,lm2,lm3))

X <- select(d,siteno,cat_area_sqkm:tot_wildfire_2011) %>%
  gather(variable,value,-siteno) %>%
  group_by(siteno) %>%
  filter(value > 0) %>% # remove all zero columns
  ungroup() %>% 
  spread(variable, value) %>%
  select_if(~!any(is.na(.))) %>% # remove NA columns
  select(-siteno) %>%
  mutate_all(funs(scale)) %>%
  as.matrix()

# create train and test set
set.seed(3)
ind <- sample(2, nrow(y), replace=T, prob=c(0.70, 0.30))

Xtrain <- X[ind==1,]
ytrain <- y[ind==1,]
Xtest <- X[ind==2,]
ytest <- y[ind==2,]

# Build model for single output
# Add layers to the model
model <- keras_model_sequential()  %>% 
  layer_dense(units = 64, input_shape=dim(Xtrain)[2], activation='relu') %>% 
  layer_dense(units = 64, activation='relu') %>% 
  layer_dense(units = 64, activation='relu') %>%
  layer_dense(units = 1, activation='relu') %>%
  compile(optimizer = "rmsprop",loss="mse",metrics="mae")
  
model_fit <- model %>% 
  fit(Xtrain, 
      ytrain[,1], 
      epochs=130, 
      batch_size = 5, 
      validation_split = 0.2, 
      callbacks = callbacks,
      verbose=0)

plot(model_fit)

# make predictions
yhat <- predict(model, Xtest)
preds <- data.frame(y=ytest[,1],yhat)
plot(preds)

# Build model for multiple outputs
input <- layer_input(shape=dim(X)[2],name="basinchars")

base_model <- input  %>%
  layer_dense(units = 64, activation='relu') %>% 
  layer_dense(units = 64, activation='relu') %>%
  layer_dense(units = 64, activation='relu') 

lm1_pred <- base_model %>% layer_dense(units = 1, name="lm1") 
lm2_pred <- base_model %>% layer_dense(units = 1, name="lm2") 
lm3_pred <- base_model %>% layer_dense(units = 1, name="lm3") 

model <- keras_model(input,list(lm1_pred,lm2_pred,lm3_pred)) %>%
  compile(optimizer = "rmsprop",
          loss="mse",
          metrics="mae",
          loss_weights=c(1,1,1))


model_fit <- model %>% 
  fit(x=X,
      y=list(y[,1],y[,2],y[,3]), 
      epochs=20, 
      batch_size = 64, 
      validation_split = 0.5, 
      verbose=0)

plot(model_fit,metrics=c("loss","lm1_loss","lm2_loss","lm3_loss"),theme_bw=T)


  
