library(tidyverse)
library(caret)
library(doMC)
library(feather)

gages <- read_feather("data/gage/all_gage_covariates2.feather") %>%
  select(acc_hdens:acc_rdx) %>%
  mutate(origin="gaged") %>%
  select(origin,everything())

huc12s <- read_feather("data/huc12/all_huc12_covariates2.feather") %>%
  filter(acc_hdens != -9999) %>%
  select(acc_hdens:acc_rdx) %>%
  select(-perennial_ice_snow) %>%
  mutate(origin="ungaged") %>%
  select(origin,everything())

# combine and ramdomly shuffle
all <- gages %>%
  bind_rows(huc12s) %>%
  mutate(origin=factor(origin)) %>%
  mutate_at(vars(-origin),funs(as.vector(scale(.)))) %>%
  select_if(~!any(is.na(.))) %>%
  sample_n(nrow(.))

# create train and test set
set.seed(0)
ind <- sample(2, nrow(all), replace=T, prob=c(0.70, 0.30))

train <- all[ind==1,]
test <- all[ind==2,]

# RF classification
registerDoMC(cores = 4)
control <- trainControl(method="repeatedcv", number=4, repeats=1)
set.seed(1)
mtry <- sqrt(ncol(train)-1)
tunegrid <- expand.grid(.mtry=mtry)
rf_fit <- train(origin~., 
                data=train, 
                method="rf", 
                tuneGrid=tunegrid, 
                trControl=control)

saveRDS(rf_fit, "data/rf_cov_shift.rds")

preds <- data.frame(obs = test$origin,
                    est = predict(rf_fit,test)) 

# confusionMatrix(data = preds$est, reference = preds$obs)

ggplot(preds) +
  geom_histogram(aes(est),stat="count") +
  facet_wrap(~obs,scale="free_y")

# mcc(preds$obs,preds$est)
# 
# mcc <- function (actual, predicted)
# {
#   # handles zero denominator and verflow error on large-ish products in denominator.
#   #
#   # actual = vector of true outcomes, 1 = Positive, 0 = Negative
#   # predicted = vector of predicted outcomes, 1 = Positive, 0 = Negative
#   # function returns MCC
# 
#   TP <- sum(actual == 1 & predicted == 1)
#   TN <- sum(actual == 0 & predicted == 0)
#   FP <- sum(actual == 0 & predicted == 1)
#   FN <- sum(actual == 1 & predicted == 0)
# 
#   sum1 <- TP+FP; sum2 <-TP+FN ; sum3 <-TN+FP ; sum4 <- TN+FN;
#   denom <- as.double(sum1)*sum2*sum3*sum4 # as.double to avoid overflow error on large products
# 
#   if (any(sum1==0, sum2==0, sum3==0, sum4==0)) {
#     denom <- 1
#   }
# 
#   mcc <- ((TP*TN)-(FP*FN)) / sqrt(denom)
#   return(mcc)
# }
