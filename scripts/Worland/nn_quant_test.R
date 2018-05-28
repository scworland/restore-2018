
d <- read_feather("data/gage/all_gage_data.feather")

Y <- select(d,area = tot_basin_area,f0.02) %>%
  mutate(f0.02 = log10((f0.02+0.001)/area)) %>%
  select(-area) %>%
  as.matrix()

X <- select(d,ppt_mean:tot_rdx) %>%
  select(-tot_basin_area) %>%
  mutate_all(funs(as.vector(scale(.)))) %>%
  select_if(~!any(is.na(.))) %>% 
  as.matrix()

build_model <- function(){
  model <- keras_model_sequential()  %>%
    layer_dense(units = 33,activation="tanh",input_shape=dim(X)[2]) %>%
    layer_dense(units = 30,activation="tanh") %>%
    layer_dense(units=1)
  
  model %>% 
    compile(optimizer = "rmsprop",
            loss="mse",
            metrics="mae")
  
  return(model)
}

k <- 5
N <- nrow(d)
indices <- sample(1:N)
folds <- cut(indices, breaks = k, labels = FALSE)
pb <- txtProgressBar(min = 0, max = k, style = 3)

mse_all <- NULL
preds_all <- NULL
for (i in 1:k) {
  
  setTxtProgressBar(pb, i)
  
  ki <- which(folds == i)
  n <- length(ki)
  
  # training set
  Xtrain <- X[-ki,]
  Ytrain <- Y[-ki,] 
  
  # validation set
  Xtest <- X[ki,]
  Ytest <- Y[ki,] 
  
  # build and fit model
  model <- build_model()
  model_fit <- model %>% 
    fit(x=Xtrain,
        y=Ytrain, 
        epochs=250, 
        batch_size = 100,
        validation_data = list(Xtest, Ytest),
        verbose=0)
  
  mse_all <- rbind(mse_all, model_fit$metrics$loss)
  
  yhat <- predict(model, Xtest) 
  
  preds <- data.frame(site_no=d$site_no[ki],
                      decade=d$decade[ki],
                      area=d$tot_basin_area[ki],
                      yhat=yhat,
                      obs=Ytest)
  
  preds_all <- rbind(preds_all,preds)
  
  close(pb)
}

average_mse <- data.frame(epoch = seq(1:ncol(mse_all))) %>%
  mutate(val_mse = sqrt(apply(mse_all, 2, mean)))

ggplot(average_mse) +
  geom_line(aes(x = epoch, y = val_mse))

est_obs <- preds_all %>%
  mutate(yhat =round((10^(yhat)*area)-0.001,3),
         obs = round((10^(obs)*area)-0.001,3)) 

ggplot(est_obs) + geom_point(aes(obs,yhat))


