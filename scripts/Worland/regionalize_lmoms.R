
library(pacman)
pacman::p_load(tidyverse,feather,lmomco,MultivariateRandomForest)
source("scripts/worland/utils.R")

d <- read_feather("data/gage/all_gage_data.feather")

Y <- select(d,L1:T6) 

X <- select(d,ppt_mean:tot_rdx) %>%
  mutate_all(funs(as.vector(scale(.)))) %>%
  select_if(~!any(is.na(.)))

# leave one-out-CV L-moments
model_data <- data.frame(Y,X)

fill <- rep(NA,nrow(d))
pred_lmoms <- data.frame(site_no=d$site_no,decade=d$decade,
                         L1=fill,L2=fill,T3=fill)

for (i in 1:nrow(d)){
  train <- model_data[-i,]
  
  l1 <- lm(log10(L1)~., select(train,L1,ppt_mean:tot_rdx))
  l2 <- lm(log10(L2)~., select(train,L2,ppt_mean:tot_rdx))
  t3 <- lm(T3~., select(train,T3,ppt_mean:tot_rdx))
  
  test <- model_data[i,]
  
  pred_lmoms$L1[i] <- 10^predict(l1,test)
  pred_lmoms$L2[i] <- 10^predict(l2,test)
  pred_lmoms$T3[i] <- predict(t3,test)
}

obs_lmoms <- select(Y,L1,L2,T3)

compare <- data.frame(ff = nonexceeds()) %>%
  mutate(est = sw_lmoms2(pred_lmoms[1,]),
         obs = sw_lmoms2(obs_lmoms[1,])) %>%
  gather(fdc,value,-ff)
                      

ggplot(compare) +
  geom_line(aes(ff,value,linetype=fdc)) +
  scale_y_continuous(trans='log2') +
  labs(x="exceedence probability",y="Q") +
  ggtitle("site #02293986") +
  theme_bw()


# old ----
d <- read_feather("data/gage/lmom3.feather")

Y <- select(d,lm1,lm2,lm3) %>%
  as.matrix()

X <- select(d,siteno,cat_area_sqkm:tot_wildfire_2011) %>%
  gather(variable,value,-siteno) %>%
  group_by(siteno) %>%
  filter(value > 0) %>% 
  ungroup() %>%
  spread(variable, value) %>%
  select_if(~!any(is.na(.))) %>%
  select(-siteno) %>%
  mutate_all(funs(scale)) %>%
  as.matrix()
  
# hold <- CrossValidation(X, Y, 5)

Yhat <- build_forest_predict(X, Y, n_tree=100, m_feature=20, min_leaf=5, X)

output <- data.frame(Y,Yhat)

