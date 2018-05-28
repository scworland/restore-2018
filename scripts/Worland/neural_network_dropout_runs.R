
library(keras)
library(tidyverse)
library(feather)
library(forcats)
library(doParallel)

setwd("~/Documents/Restore")

source("scripts/Worland/clean/utils.R")

# setwd("~/Documents/Restore")
d <- read_feather("data/gage/all_gage_data2.feather") 

Y <- select(d,f0.02:f99.98) %>%
  mutate_all(funs(.+2)) %>%
  mutate_all(funs(log10))

X <- select(d,lon,lat,acc_hdens:statsgo) %>%
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

cl <- makeCluster(5) # change
registerDoParallel(cl)
run <- foreach(i=1:100) %dopar% {
  cv_results <- k_foldcv(k=10,epochs=200,batch_size=30,X=X,Y=Y,data=d,pls=FALSE)
  cv_results$est_obs
}
stopCluster(cl)

saveRDS(run, "data/sandbox/dropout_20.rds")

run <- readRDS("data/sandbox/dropout_100.rds")

sites_decade <- read_feather("data/sandbox/combined_est.feather") %>%
  distinct(site_no,decade)

runs_combined <- run %>%
  purrr::reduce(left_join, by = c("site_no","decade","variable","obs")) %>%
  select(site_no,decade,variable,obs,everything()) %>%
  setNames(c("site_no","decade","variable","obs",paste0("run_",1:length(run)))) %>%
  gather(run,value,-site_no,-decade,-variable,-obs) %>%
  mutate(obs = (10^obs)-2,
         value = (10^value)-2,
         value = ifelse(value<0,0,value),
         variable = as.numeric(substring(variable, 2))) %>%
  group_by(site_no,decade,variable) %>%
  summarize(obs = mean(obs),
            mu = mean(value),
            min = min(value),
            max = max(value)) %>%
  right_join(sites_decade,by=c("site_no","decade"))

ggplot(runs_combined) +
  geom_line(aes(variable,obs)) +
  geom_line(aes(variable,min),linetype="dashed") +
  geom_line(aes(variable,max),linetype="dashed") +
  geom_line(aes(variable,mu),linetype="dotted") +
  facet_wrap(~site_no,scales="free_y") +
  labs(x="non-exceedance probability",y="Q") +
  theme_bw() +
  scale_y_log10()



