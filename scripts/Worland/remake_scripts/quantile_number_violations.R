library(keras)
library(tidyverse)
library(feather)
library(broom)
library(zoo)
source("scripts/Worland/utils.R")
source("scripts/Worland/keras_funs.R")

# setwd("~/Documents/Restore")
d <- remake::fetch('gage_all')

Y <- select(d,f0.02:f99.98) %>%
  mutate_all(funs(.+2)) %>%
  mutate_all(funs(log10))

X <- select(d,lon,lat,acc_hdens:statsgo,-perennial_ice_snow,-area_sqkm) %>%
  mutate_at(vars(bedperm:statsgo),funs(as.numeric(as.factor(.)))) %>%
  mutate_all(funs(as.vector(scale(.)))) %>%
  select_if(~!any(is.na(.))) %>% 
  as.matrix()

f1 <- c("f0.02","f0.05","f0.1","f0.2","f0.5","f01","f02","f05",
        "f10","f20","f25","f30","f40","f50","f60","f70","f75",
        "f80","f90","f95","f98","f99","f99.5","f99.8","f99.9",
        "f99.95","f99.98")

f2 <- c("f0.02","f0.1","f0.5","f02","f05","f10","f20","f25",
        "f30","f40","f50","f60","f70","f75","f80","f90","f95",
        "f98","f99.5","f99.9","f99.98")

f3 <- c("f0.02","f0.5","f05","f10","f20", "f30","f40","f50",
        "f60","f70","f80","f90","f95","f99.5","f99.98")

f4 <- c("f0.02","f0.5","f05","f10","f20", "f30","f40","f50",
        "f60","f70","f80","f90","f95","f99.5","f99.98")

f5 <- c("f05","f10","f20", "f30","f40","f50","f60","f70",
        "f80","f90","f95")

f6 <- c("f05","f20", "f30","f40","f50","f60","f70","f80","f95")

f_all <- list(f1,f2,f3,f4,f5,f6)

viol_all <- NULL
for (i in 1:length(f_all)) {
  
  print(paste0("Analyzing Scenario ",i," out of ", length(f_all)))
  
  # subset Y matrix
  Y2 <- Y[,f_all[[i]]]
  
  # Build model function for multiple outputs
  build_model <- function(){
    input <- layer_input(shape=dim(X)[2],name="basinchars")
    
    base_model <- input  %>%
      layer_dense(units = 40,activation="relu") %>%
      layer_dropout(rate=0.1) %>%
      layer_dense(units = 20,activation="relu") %>%
      layer_dropout(rate=0.1) %>%
      layer_dense(units = 20,activation="relu") %>%
      layer_dropout(rate=0.1) 
    
    for(j in 1:dim(Y2)[2]){
      y <- colnames(Y2)[j]
      outstring <- paste0(
        sprintf("%s <- base_model %%>%%", y), 
        sprintf(" layer_dense(units = 1, activation='relu', name='%s')",y)
      )
      eval(parse(text=outstring))
    }
    
    Ylist <- paste0("list(",paste(colnames(Y2),sep="",collapse=","),")")
    model <- keras_model(input,eval(parse(text=Ylist))) %>%
      compile(optimizer = "rmsprop",
              loss="mse",
              metrics="mae")
    
    return(model)
  }
  
  # K-fold cross validation
  epochs <- 200
  cv_results <- k_foldcv(k=5,epochs=epochs,batch_size=100,Y=Y2,X=X,data=d,pls=FALSE)
  
  viol_count <- data.frame(t(cv_results$violations)) %>%
    set_names(paste0("kfold",1:ncol(.))) %>%
    mutate(epoch=1:epochs) %>%
    gather(fold,value,-epoch) %>%
    mutate(quants_est = paste0(ncol(Y2), "_quantiles"))
  
  viol_all <- rbind(viol_all,viol_count)
}

viol_count <- viol_all %>%
  group_by(epoch,quants_est,quants_est) %>%
  summarize(min=min(value),
            max=max(value),
            mu=mean(value))


ggplot(viol_count) +
  geom_line(aes(x=epoch,y=mu,color=quants_est, group = quants_est)) +
  #scale_linetype_manual(values = rep("solid",5),guide = FALSE) +
  coord_cartesian(xlim=c(0,150)) +
  #scale_y_log10() +
  theme_bw() +
  labs(color="#quantiles est", fill = "#quantiles est", y="number of violations")

