library(tidyverse)
library(keras)
library(lime)
library(corrr)
library(forcats)
library(tidyquant)

# number of observations
N <- 1000

# synthetic multiple outputs
y <- rnorm(N,20,5) # smallest

# synthetic nonlinear predictors
x1 <- (sqrt(y) + rnorm(N,0,1))^2
x2 <- log10((x1 + y)^3)

X <- data.frame(x1,x2) %>% # combine
  mutate_all(funs(scale)) %>% # scale
  as.matrix() # to matrix

# create train and test set
set.seed(0)
ind <- sample(2, N, replace=T, prob=c(0.50, 0.50))

Xtrain <- X[ind==1,]
Ytrain <- y[ind==1]

Xtest <- X[ind==2,]
Ytest <- y[ind==2] 

model <- keras_model_sequential() %>%
  layer_dense(units = 2, input_shape = dim(X)[2]) %>%
  layer_dense(units = 30,activation="relu") %>%
  layer_dropout(rate = 0.1) %>%
  layer_dense(units = 1) %>%
  compile(optimizer = "rmsprop",
          loss="mse",
          metrics="mae")

model_fit <- model %>% 
  fit(x=Xtrain,
      y=Ytrain, 
      epochs=100, 
      batch_size = 25, 
      validation_split = 0.2,
      verbose=0)

Xtrain <- as.tibble(Xtrain)
Xtest <- as.tibble(Xtest)

# setup lime::predict_model() function for keras
predict_model.keras.models.Sequential <- function(x,newdata,type...) {
  pred <- predict(x, x=as.matrix(newdata)) %>%
    data.frame()
}

# setup lime::predict_model() function for keras
# predict_model.keras.engine.training.Model <- function(x,newdata,type...) {
#   pred <- predict(x, x=as.matrix(newdata)) %>%
#     data.frame()
# }

# test predict model function
predict_model(x = model, newdata = Xtest, type='raw') %>%
  tibble::as.tibble()

# run lime() on training set
explainer <- lime::lime(x = Xtrain, 
                        model = model, 
                        bin_continuous = FALSE)

# run explain() on the explainer on first 5 rows
explanation <- lime::explain(x = Xtest[1:50,], 
                             explainer = explainer, 
                             n_features = 4,
                             kernel_width = 0.5)

summary_explain <- explanation %>%
  group_by(feature_desc) %>%
  summarize(low=quantile(feature_weight,0.1),
            high=quantile(feature_weight,0.9),
            mu=mean(feature_weight))

# lime::plot() options
plot_features(explanation)
plot_explanations(explanation)

# correlation analysis
corrr_analysis <- Xtrain %>%
  mutate(y = Ytrain) %>%
  correlate() %>%
  focus(y) %>%
  rename(feature = rowname) %>%
  arrange(abs(y)) %>%
  mutate(feature = as.factor(feature))

corrr_analysis %>%
  ggplot(aes(x = y, y = fct_reorder(feature, desc(y)))) +
  geom_point() +
  # Positive Correlations - Contribute to churn
  geom_segment(aes(xend = 0, yend = feature), 
               color = palette_light()[[2]], 
               data = corrr_analysis %>% filter(y > 0)) +
  geom_point(color = palette_light()[[2]], 
             data = corrr_analysis %>% filter(y > 0)) +
  # Negative Correlations - Prevent churn
  geom_segment(aes(xend = 0, yend = feature), 
               color = palette_light()[[1]], 
               data = corrr_analysis %>% filter(y < 0)) +
  geom_point(color = palette_light()[[1]], 
             data = corrr_analysis %>% filter(y < 0)) +
  theme_bw()
