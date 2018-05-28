
library(keras)
library(tensorflow)
library(tidyverse)
library(dplyr)

# number of observations
N <- 1000

# synthetic multiple outputs
y1 <- rnorm(N,20,5) # smallest
y2 <- y1 + rnorm(N,7,2) # middle
y3 <- y2 + rnorm(N,7,2) # largest
Y <- data.frame(y1,y2,y3) # combine
Y <- Y * rexp(N,2) # randomly scale rows

# synthetic nonlinear predictors
x1 <- (sqrt(y1) + rnorm(N,0,1))^2
x2 <- log10((x1 + y2)^3)
x3 <- sqrt(abs(y3-y2))
x4 <- sqrt((x1 + y2)^2)

X <- data.frame(x1,x2,x3,x4) %>% # combine
  mutate_all(funs(scale)) %>% # scale
  as.matrix() # to matrix

# pairs(cbind(Y,X))

# pls_fit <- plsreg2(X,Y,comps=2)

# create train and test set
set.seed(0)
ind <- sample(2, N, replace=T, prob=c(0.50, 0.50))

Xtrain <- X[ind==1,]
Ytrain <- Y[ind==1,] %>%
  as.matrix() %>%
  split(.,col(.)) %>%
  unname() %>%
  keras_array()

Xtest <- X[ind==2,]
Ytest <- Y[ind==2,] %>%
  as.matrix() %>%
  split(.,col(.)) %>%
  unname() %>%
  keras_array()


# add covariates
input <- layer_input(shape=dim(X)[2],name="covars")

# add hidden layers
base_model <- input  %>%
  layer_dense(units = 4, activation='relu') 

# add output 1 
yhat1 <- base_model %>% 
  layer_dense(units = 1, name="yhat1") 

# add output 2
yhat2 <- base_model %>% 
  layer_dense(units = 1, name="yhat2") 

# add output 3
yhat3 <- base_model %>% 
  layer_dense(units = 1, name="yhat3")

# caclulate squared weights for mse loss
weight_fun <- function(x){mean(Y$y1/x)}
weights <- as.numeric(apply(Y,2,weight_fun))^2

# define custom callback class
violationHistory <- R6::R6Class(
  "violationHistory",
  inherit = KerasCallback,
  
  public = list(
    val = NULL,
    
    on_epoch_end = function(epoch, logs = list()) {
      preds <- do.call(cbind, model$predict(Xtest))
      value <-
        sum(apply(preds, 1, function(x) {
          sum(cummax(x) != x)
        }))
      self$val <- c(self$val, value)
    }
  )
)

# custom loss function
custom_loss <- function(y_true, y_pred){
  preds <- do.call(cbind, model$predict(Xtest))
  assign("val", violation_count$val, envir = .GlobalEnv)
  mse <- k_mean(k_square(y_true-y_pred)) 
  val <- k_cast(val,'float32')
  return(val+mse)
}

# build multi-output model
violation_count <- violationHistory$new()
model <- keras_model(input,list(yhat1,yhat2,yhat3)) %>%
  compile(optimizer = "rmsprop",
          loss='mse',
          metrics='mae',
          loss_weights=weights)

# fit model (intentionally underfit)
model_fit <- model %>% 
  fit(x=Xtrain,
      y=Ytrain, 
      epochs=100, 
      batch_size = 25, 
      validation_data = list(Xtest,Ytest),
      verbose=0,
      callbacks=violation_count)

data.frame(epoch=1:100,violations=violation_count$val) %>%
  ggplot(aes(epoch,violations)) +
  geom_line() +
  #geom_point() +
  theme_bw()

# plot overall loss
data.frame(loss=model_fit$metrics$loss,
           epoch=1:100) %>%
  ggplot() + 
  geom_line(aes(epoch,loss)) +
  theme_bw()

# predict values for test set
Yhat <- predict(model, Xtest) %>%
  data.frame() %>%
  setNames(colnames(Y)) 

predictions <- model %>% 
  predict(as.matrix(Xtest)) %>% 
  data.frame() %>%
  setNames(colnames(Y)) 

# plot est vs obs
Yhat %>%
  gather(variable,yhat) %>%
  mutate(obs = gather(Ytest,var,obs)$obs) %>%
  ggplot() +
  geom_point(aes(obs,yhat)) +
  facet_wrap(~variable,ncol=3) +
  geom_abline(linetype="dashed") +
  theme_bw()

# sum monotonic violations 
sum(apply(Yhat,1,function(x){sum(cummax(x)!=x)}))
sum(apply(Ytest,1,function(x){sum(cummax(x)!=x)}))





# custom loss function
# loss_custom <- function(y_true, y_pred){
#   val <- k_placeholder(dtype = 'float32')
#   mse <- k_placeholder(dtype = 'float32')
#   n <- k_placeholder(dtype = 'float32')
#   k <- k_placeholder(dtype = 'float32')
#   n <- k_shape(y_pred)[[1]]
#   k <- k_shape(y_pred)[[2]]
#   y_pred <- array_reshape(y_pred,dim=c(k,n))
#   y_true <- array_reshape(y_true,dim=c(k,n))
#   val <- k_sum(k_cast(k_less(y_pred[2:k]-y_pred[1:k-1],0),'float32'))
#   mse <- k_mean(k_square(y_true-y_pred)) 
#   return(mse + val)
# }
