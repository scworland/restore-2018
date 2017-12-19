
library(tidyverse)
library(lmomco)
library(lubridate)
library(feather)
source("scripts/worland/utils.R")

gt <- read_feather("data/gage/gage_time.feather")
gs <- read_feather("data/gage/gage_static.feather")

n = length(unique(gt$siteno))

lmom3 <- gt %>%
  select(siteno,date,Q) %>%
  mutate(Q=log(Q+0.001)) %>%
  group_by(siteno) %>% 
  do(lms=lmoms(.$Q, nmom=3)$lambdas) %>%
  ungroup() %>%
  unnest(lms) %>%
  mutate(names = rep(c("lm1","lm2","lm3"),n)) %>%
  spread(names,lms) %>%
  left_join(gs,by="siteno")

write_feather(lmom3,"data/gage/lmom3.feather")

test <- gt %>%
    select(siteno,date,Q) %>%
    filter(siteno %in% unique(siteno)[1:3]) %>%
    group_by(siteno) %>%    
    mutate(p = sw_efdc(log(Q)),
         p2 = sw_lmoms(log(Q),type="pe3"))

ggplot(test) + 
  geom_line(aes(p,Q)) + 
  geom_line(aes(p2,Q),linetype="dashed",color="dodgerblue") +
  facet_wrap(~siteno, scales="free_y") +
  scale_y_continuous(trans='log2') +
  theme_bw()


# prob x < value 
prob <- pnorm(q=c(28,100), mean=50, sd=20)

# quantile associated with prob
quant <- qnorm(p=0.8,mean=50,sd=20)

x <- 5
y <- 3

alpha <- rnorm(10,2,0.5)
beta <- (y-alpha)/x

ggplot() + 
  geom_point(aes(x,y)) +
  coord_cartesian(xlim=c(0,10),ylim=c(0,10)) +
  theme_bw() +
  geom_abline(intercept=alpha,slope=beta)


