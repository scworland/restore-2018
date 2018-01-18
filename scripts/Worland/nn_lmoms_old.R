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

# ------------------------
# K-fold cross validation
k <- 10
indices <- sample(1:nrow(d))
folds <- cut(indices, breaks = k, labels = FALSE)

mse_all <- NULL
obs_all <- NULL
yhat_all <- NULL
for (i in 1:k) {
  
  if(i %% 10==0) {
    # Print every 10th iteration
    print(paste0("iteration: ", i))
  }
  
  ki <- which(folds == i)
  n <- length(ki)
  
  # training set
  Xtrain <- X[-ki,]
  Ytrain <- keras_array(list(Y[-ki,1],Y[-ki,2],Y[-ki,3]))
  
  # validation set
  Xtest <- array_reshape(X[ki,],c(n,34))
  Ytest <- keras_array(list(Y[ki,1],Y[ki,2],Y[ki,3]))
  
  # build and fit model
  model <- build_model()
  model_fit <- model %>% 
    fit(x=Xtrain,
        y=Ytrain, 
        epochs=40, 
        batch_size = 250,
        validation_data = list(Xtest, Ytest), 
        verbose=0)
  
  mse_all <- rbind(mse_all, model_fit$metrics$val_loss)
  obs_all <- rbind(obs_all,data.frame(Y[ki,]))
  
  yhat <- predict(model, Xtest) %>%
    data.frame() %>%
    setNames(c("L1","L2","T3")) 
  
  yhat_all <- rbind(yhat_all,yhat)
}

average_mse <- data.frame(
  epoch = seq(1:ncol(mse_all)),
  validation_mse = apply(mse_all, 2, mean))

ggplot(average_mse, aes(x = epoch, y = validation_mse)) + geom_line()
