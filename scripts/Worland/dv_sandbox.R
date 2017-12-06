
library(tidyverse)
library(dataRetrieval)
library(lmomco)
library(lubridate)
library(feather)
source("scripts/worland/utils.R")

gt <- read_feather("data/gage/gage_time.feather")
gs <- read_feather("data/gage/gage_static.feather")

test <- gt %>%
  select(siteno,date,Q) %>%
  filter(siteno %in% unique(siteno)[10:21]) %>%
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


