library(tidyverse)
library(feather)
library(keras)
library(lime)
library(corrr)

source("scripts/Worland/utils.R")

# load data
d <- remake::fetch('gage_all')

Y <- select(d,f0.02:f99.98) %>%
  mutate_all(funs(.+2)) %>%
  mutate_all(funs(log10)) %>%
  as.matrix()

X <- select(d,lon,lat,acc_hdens:statsgo,-perennial_ice_snow,-area_sqkm) %>%
  mutate_at(vars(bedperm:statsgo),funs(as.numeric(as.factor(.)))) %>%
  mutate_all(funs(as.vector(scale(.)))) %>%
  select_if(~!any(is.na(.))) %>% 
  as.matrix()

# create train and test set
set.seed(0)
ind <- sample(2, nrow(d), replace=T, prob=c(0.50, 0.50))

Xtrain <- X[ind==1,]
Ytrain <- Y[ind==1,]

Xtest <- X[ind==2,]
Ytest <- Y[ind==2,] 

explanation_all <- NULL
for(i in 1:ncol(Y)){
  
  print(paste0("Calculating LIME for quantile #",i," out of ", ncol(Y)))
  
  model <- keras_model_sequential() %>%
    layer_dense(units = dim(X)[2], input_shape = dim(X)[2]) %>%
    layer_dense(units = 30,activation="relu") %>%
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 1) %>%
    compile(optimizer = "rmsprop",
            loss="mse",
            metrics="mae")
  
  model_fit <- model %>% 
    fit(x=Xtrain,
        y=Ytrain[,i], 
        epochs=100, 
        batch_size = 25, 
        validation_split = 0.2,
        verbose=0)
  
  # setup lime::predict_model() function for keras
  predict_model.keras.models.Sequential <- function(x,newdata,type...) {
    pred <- predict(x, x=as.matrix(newdata)) %>%
      data.frame()
  }
  
  # run lime() on training set
  explainer <- lime::lime(x = as.tibble(Xtrain), 
                          model = model, 
                          bin_continuous = FALSE)
  
  # run explain() on the explainer on first 5 rows
  explanation <- lime::explain(x = as.tibble(Xtest[1:100,]), 
                               explainer = explainer, 
                               n_features = 10,
                               feature_select="highest_weights",
                               kernel_width = 0.5) %>%
    mutate(quantile = colnames(Y)[i])
  
  explanation_all <- rbind(explanation_all,explanation)
  
  rm(model)
}

plot_final <- explanation_all %>%
  add_count(feature) %>%
  mutate(quantile = substring(quantile, 2)) %>%
  #         as.numeric(.)/100) %>%
  filter(n > 1000)
  
ggplot(plot_final) +
  geom_boxplot(aes(x=quantile, y=feature_weight),outlier.shape = NA) +
  # geom_linerange(aes(x=quantile,ymin=low,ymax=high)) +
  # geom_line(aes(x=quantile,y=mu)) +
  facet_wrap(~feature, scales="free_y") + 
  geom_hline(yintercept=0,linetype="dashed") +
  coord_cartesian(ylim=c(-0.85,0.6)) +
  theme_bw()
