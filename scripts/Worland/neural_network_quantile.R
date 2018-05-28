
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

X <- select(d,lon,lat,acc_hdens:statsgo,-perennial_ice_snow,-area_sqkm) %>%
  mutate_at(vars(bedperm:statsgo),funs(as.numeric(as.factor(.)))) %>%
  mutate_all(funs(as.vector(scale(.)))) %>%
  select_if(~!any(is.na(.))) %>% 
  as.matrix()

# Build model function for multiple outputs
build_model <- function(){
  input <- layer_input(shape=dim(X)[2],name="basinchars")
  
  base_model <- input  %>%
    layer_dense(units = 40,activation="relu") %>%
    layer_dropout(rate=0.1) %>%
    layer_dense(units = 30,activation="relu") %>%
    layer_dropout(rate=0.1) 
  
  for(i in 1:dim(Y)[2]){
    y <- colnames(Y)[i]
    outstring <- paste0(
      sprintf("%s <- base_model %%>%%", y), 
      sprintf(" layer_dense(units = 1, activation='relu', name='%s')",y)
    )
    eval(parse(text=outstring))
  }
  
  Ylist <- paste0("list(",paste(colnames(Y),sep="",collapse=","),")")
  model <- keras_model(input,eval(parse(text=Ylist))) %>%
    compile(optimizer = "rmsprop",
            loss="mse",
            metrics="mae")
  
  return(model)
}

# K-fold cross validation
source("scripts/Worland/keras_funs.R")

cv_results <- k_foldcv(k=2,epochs=100,batch_size=80,Y=Y,X=X,data=d,pls=FALSE)

est_obs <- cv_results$est_obs %>%
  mutate(obs = (10^obs)-2,
         nnet = round((10^nnet)-2,2),
         variable = as.numeric(substring(variable, 2)))

write_feather(est_obs,"data/gage/direct_estimates.feather")

# tail(cv_results$average_mse)
# 
# ggplot(cv_results$average_mse) +
#   geom_line(aes(x = epoch, y = val_mse))
# 
# # fit metrics
# metrics <- cv_results$est_obs %>%
#   mutate(obs = (10^obs)-2,
#          nnet = round((10^nnet)-2,2),
#          variable = as.numeric(substring(variable, 2))) %>%
#   group_by(variable,decade) %>%
#   summarize(corr = cor(obs,nnet),
#             rmse = sqrt(mean((obs-nnet)^2))/mean(obs),
#             mpe = median(abs((nnet[obs>0]-obs[obs>0])/obs[obs>0]))*100) %>%
#   gather(metric,value,-variable,-decade)
#   
# ggplot(metrics) +
#   geom_line(aes(variable,value)) +
#   facet_grid(metric~decade,scales="free_y") +
#   labs(x = "exceedance probability") +
#   theme_bw()
# 
# # plots estimated vs observed
# est_obs <- cv_results$est_obs %>%
#   mutate(obs = (10^obs)-2,
#          nnet = round((10^nnet)-2,2),
#          variable = as.numeric(substring(variable, 2))) %>%
#   group_by(site_no,decade) %>%
#   mutate(nnet = sort(nnet)) %>%
#   # mutate(nnet = ifelse(nnet<0,0,nnet) %>%
#   #          ifelse(.==cummax(.),.,NA) %>%
#   #          na.approx(.,rule=2)) %>%
#          #loess_fit = predict(loess(nnet~variable,span=0.2))) %>%
#   gather(model,value,-site_no,-decade,-variable) %>%
#   group_by(site_no) %>%
#   sample_n_groups(6)
# 
# ggplot(est_obs) +
#   geom_line(aes(variable,value,linetype=model)) +
#   scale_linetype_manual(values=c("dashed","solid","dotted")) +
#   facet_wrap(site_no~decade,scales="free_y") +
#   labs(x="EP",y="Q") +
#   theme_bw() +
#   scale_y_log10()
# 
# # constrained GAM
# # hold <- cv_results$est_obs %>%
# #   rename(nn=model) %>%
# #   left_join(scales,by=c("variable","decade")) %>%
# #   mutate(nn = (nn*sigma)+mu-2) %>%
# #   select(-mu,-sigma) %>%
# #   filter(site_no == "02352500")  %>%
# #   mutate(variable = as.numeric(substring(variable, 2))) %>%
# #   group_by(decade,site_no) %>%
# #   do(augment(scam(nn ~ s(variable, k=14,bs="mpi"), data =.))) 
# 
# # Y <- select(d,decade,f0.02:f99.98) %>%
# #   mutate_all(funs(replace(., .==0, 0.001))) %>%
# #   mutate_all(funs(log10)) 
# 
# 
hold <- data.frame(t(cv_results$violations)) %>%
  set_names(paste0("kfold",1:ncol(.))) %>%
  mutate(epoch=1:100) %>%
  gather(fold,value,-epoch) %>%
  mutate(paste0(ncol(Y2), "_quantiles"))

ggplot(hold) +
  geom_line(aes(epoch,value,group=fold),alpha=0.5) +
  theme_bw() +
  scale_y_log10()

