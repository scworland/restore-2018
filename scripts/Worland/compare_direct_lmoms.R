
library(lmomco)
library(tidyverse)
library(feather)
library(parallel)
library(gridExtra)
library(sf)
source("scripts/Worland/utils.R")

d <- read_feather("data/gage/all_gage_data2.feather") 

direct_est <- read_feather("data/gage/direct_estimates.feather") %>%
  rename(multi_nnet = nnet)

direct_est_ind <- read_feather("data/gage/direct_indv_estimates.feather") %>%
  rename(ind_nnet = nnet)

lmom_est <- read_feather("data/gage/lmom_estimates.feather")

# count number of NAs
# countNA <- aep_all %>% 
#   group_by(site_no,decade) %>% 
#   summarize(n=sum(is.na(aep))) %>% 
#   filter(n>0)

# plot FDCs
combined_est <- direct_est %>%
  left_join(lmom_est,by=c("site_no","decade","variable")) %>%
  na.omit() %>%
  gather(model,value,-site_no,-decade,-variable) %>%
  group_by(site_no,decade) %>%
  sample_n_groups(1) 


# ggplot(filter(combined_est,model=="obs")) +
#   geom_line(aes(variable,value)) +
#   geom_point(aes(variable,value),shape=21,fill="white") +
#   labs(x="non-exceedance probability",y="Q") +
#   theme_bw() +
#   scale_y_log10() 

ggplot(combined_est) +
  geom_line(aes(variable,value,color=model,linetype=model)) +
  scale_color_manual(values=c("black","slateblue3","black"),
                     labels=c("AEP","Multioutput NN","Observed"),
                     name="Estimate") +
  scale_linetype_manual(values=c("solid","solid","dashed"),
                        labels=c("AEP","Multioutput NN","Observed"),
                        name="Estimate") +
  facet_wrap(~site_no,scales="free_y") +
  labs(x="non-exceedance probability",y="Q") +
  theme_bw() +
  scale_y_log10() +
  theme(legend.position = "top")

# plot percent violations and single FDC example
p1 <- direct_est %>%
  left_join(direct_est_ind,by=c("site_no","decade","variable","obs")) %>%
  group_by(site_no,decade) %>%
  summarize(ind_violations = sum(cummax(ind_nnet)!=ind_nnet),
            multi_violations = sum(cummax(multi_nnet)!=multi_nnet)) %>%
  ungroup() %>%
  group_by(decade) %>%
  #mutate(n = length(unique(site_no))) %>%
  summarize(n = length(unique(site_no)),
            ind_violations = sum(ind_violations)/(26*n),
            multi_violations = sum(multi_violations)/(26*n)) %>%
  select(-n) %>%
  gather(variable,value,-decade) %>%
  ggplot() +
  geom_line(aes(decade,value,color=variable,group=variable)) +
  scale_color_manual(values=c("black","slateblue3"),
                     labels=c("Individual quantile NN","Multioutput quantile NN"),
                     name="Model type",
                     guide=FALSE) +
  labs(y="Percent monotonic violations") +
  scale_y_continuous(labels = scales::percent) +
  #geom_hline(yintercept=0.053,linetype="dashed") +
  theme_bw() 

p2 <- direct_est %>%
  left_join(direct_est_ind,by=c("site_no","decade","variable","obs")) %>%
  filter(site_no=="02315500" & decade=="1950") %>%
  gather(model,value,-site_no,-decade,-variable) %>%
  ggplot() +
  geom_line(aes(variable,value,color=model,linetype=model)) +
  scale_color_manual(values=c("black","slateblue3","black"),
                     labels=c("Individual NN","Multioutput NN","Observed"),
                     name="Estimate") +
  scale_linetype_manual(values=c("solid","solid","dashed"),
                        labels=c("Individual NN","Multioutput NN","Observed"),
                        name="Estimate") +
  labs(x="non-exceedance probability",y="Q") +
  theme_bw() +
  scale_y_log10() +
  theme(legend.position = c(0.8,0.25)) 

grid.arrange(p1,p2,nrow=1)

# plot metrics
A = function(yhat,y){2*sd(yhat)*sd(y)*(1-cor(yhat,y))}
B = function(yhat,y){(sd(yhat) - sd(y))^2}
C = function(yhat,y){(mean(yhat) - mean(y))^2}
MSE = function(yhat,y){A(yhat,y) + B(yhat,y) + C(yhat,y)}

metrics <- direct_est %>%
  left_join(direct_est_ind,by=c("site_no","decade","variable","obs")) %>%
  #mutate(aep=round(aep,3)) %>%
  gather(model,value,-site_no,-decade,-variable,-obs) %>%
  group_by(variable,decade,model) %>%
  summarize(var = B(value,obs)/mean(obs),
            bias = C(value,obs)/mean(obs),
            corr = A(value,obs)/mean(obs)) %>%
  gather(metric,value,-variable,-decade,-model)

metrics <- direct_est %>%
  left_join(lmom_est,by=c("site_no","decade","variable")) %>%
  filter(multi_nnet > 0 & obs > 0 & aep > 0) %>%
  gather(model,value,-site_no,-decade,-variable,-obs) %>%
  mutate(rel_error = round((value-obs)/obs,2),
         error = ifelse(rel_error>1,1,rel_error),
         model = ifelse(model=="aep","AEP","Multioutput NN"),
         variable = as.factor(variable))

ggplot(metrics) + 
  geom_boxplot(aes(variable,error),fill="grey70",outlier.shape = NA,alpha=0.5) +
  coord_cartesian(ylim=c(-1.1,1.1)) +
  geom_hline(yintercept=0,linetype="dashed") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~model,nrow=1) +
  labs(x="nonexceedance probability",y="relative error") +
  theme_bw()
  
  
metrics <- direct_est %>%
  left_join(lmom_est,by=c("site_no","decade","variable")) %>%
  #mutate(aep=round(aep,3)) %>%
  gather(model,value,-site_no,-decade,-variable,-obs) %>%
  group_by(variable,decade,model) %>%
  summarize(corr = cor(obs,value,use = "pairwise.complete.obs"),
            rmse = sqrt(mean((obs-value)^2,na.rm=T))/mean(obs),
            mpe = median(abs((value[obs>0]-obs[obs>0])/obs[obs>0]),na.rm=T)*100) %>%
  gather(metric,value,-variable,-decade,-model)

ggplot(metrics) +
  geom_line(aes(variable,value,linetype=model)) +
  scale_linetype_manual(values=c("longdash","solid"),labels=c("AEP","Multioutput NN")) +
  facet_grid(metric~decade,scales="free_y") +
  labs(x = "exceedance probability") +
  #scale_y_log10() +
  theme_bw() +
  theme(legend.position = "top")

# map errors
boundary <- read_sf("data/shapefiles/simple_bounds/simple_bounds.shp")

direct_error <- direct_est %>%
  left_join(lmom_est,by=c("site_no","decade","variable")) %>%
  filter(variable==50) %>%
  filter(multi_nnet > 0 & obs > 0 & aep > 0) %>%
  gather(model,value,-site_no,-decade,-variable,-obs) %>%
  mutate(rel_error = round((value-obs)/obs,2),
         model = ifelse(model=="aep","AEP","Multioutput NN")) %>%
  group_by(site_no,model) %>%
  summarize(error = mean(rel_error)) %>%
  left_join(select(d,site_no,lon,lat),by="site_no") %>%
  distinct(site_no,model,.keep_all=TRUE) %>%
  mutate(error = ifelse(error>1,1,error))

regions = c("texas","alabama","florida","mississippi",
            "louisiana","georgia","arkansas","tennessee",
            "oklahoma")

states <- subset(map_data("state"),region %in% regions)

ggplot() +
  geom_polygon(data=states, aes(x=long, y=lat, group=group),fill="white",color="grey45") +
  geom_sf(data=boundary,fill="grey",alpha=0.8) +
  coord_sf(crs=st_crs(boundary)) +
  geom_point(data=direct_error, aes(x=lon, y=lat,color=error),
             size=2, alpha=1, shape=4) +
  scale_color_gradient2(low="blue",high="red",mid="white") +
  facet_wrap(~model,nrow=2) +
  #facet_wrap(~variable) +
  # coord_equal() +
  theme_bw()

# sandbox
hold <- direct_est %>%
  group_by(site_no,decade) %>%
  mutate(multi_nnet = ifelse(multi_nnet < 0, 0, multi_nnet),
         viol = ifelse(cummax(multi_nnet)!=multi_nnet,1,0)) %>%
  arrange(site_no,decade,variable) %>%
  filter(viol==1)

hold2 <- direct_est %>%
  filter(site_no %in% hold$site_no & decade %in% hold$decade) %>%
  group_by(site_no,decade) %>%
  mutate(multi_nnet = ifelse(multi_nnet < 0, 0, multi_nnet),
         viol = ifelse(cummax(multi_nnet)!=multi_nnet,1,0),
         dif = c(0,diff(multi_nnet))) %>%
  arrange(site_no,decade,variable)

hold3 <- direct_est %>%
  group_by(decade) %>%
  summarize(n = length(unique(site_no)))

