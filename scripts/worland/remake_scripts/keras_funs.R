

sw_k_foldcv <- function(build_model,k=10,epochs=30,batch_size=250,X=X,Y=Y,data=d){
  
  #start the clock and progress bar
  ptm <- proc.time()
  pb <- txtProgressBar(min = 0, max = k, style = 3)
  
  k <- k
  N <- nrow(data)
  set.seed(1)
  indices <- sample(1:N)
  folds <- cut(indices, breaks = k, labels = FALSE)
  
  mse_all <- NULL
  obs_all <- NULL
  yhat_all <- NULL
  null_all <- NULL
  site_no_all <- NULL
  pls_all <- NULL
  violation_all <- NULL
  for (i in 1:k) {
    
    if(i %% 2==0) {
      # print progress every 2nd iteration
      setTxtProgressBar(pb, i)
    }
    
    ki <- which(folds == i)
    n <- length(ki)

    # training set
    Xtrain <- X[-ki,]
    Ytrain <- Y %>% 
      slice(-ki) %>%
      as.matrix() %>%
      split(.,col(.)) %>%
      unname() %>%
      keras_array()
    
    # validation set
    Xtest <- array_reshape(X[ki,],c(n,ncol(X)))
    Ytest <- Y %>% 
      slice(ki) %>%
      as.matrix() %>%
      split(.,col(.)) %>%
      unname() %>%
      keras_array()
    
    # build and fit model
    model <- build_model()
    
    # define custom callback class
    violationHistory <- R6::R6Class(
      "violationHistory",
      inherit = KerasCallback,
      
      public = list(
        val = NULL,
        
        on_epoch_end = function(epoch, logs = list()) {
          preds <- do.call(cbind, model$predict(Xtest))
          value <-
            sum(apply(preds, 1, function(x) {
              sum(cummax(x) != x)
            }))
          self$val <- c(self$val, value)
        }
      )
    )

    # fit model
    if(length(Ytest)==1){
      model_fit <- model %>% 
        fit(x=Xtrain,
            y=Ytrain, 
            epochs=epochs, 
            batch_size = batch_size,
            validation_data = list(Xtest, Ytest), 
            verbose=0)
    }else{
    violation_count <- violationHistory$new()
    model_fit <- model %>% 
      fit(x=Xtrain,
          y=Ytrain, 
          epochs=epochs, 
          batch_size = batch_size,
          validation_data = list(Xtest, Ytest), 
          verbose=0)
    
    #violation_all <- rbind(violation_all,violation_count$val)
    }
    
    mse_all <- rbind(mse_all, model_fit$metrics$val_loss)
    
    if(k < N){
      obs_all <- rbind(obs_all,data.frame(Y[ki,]))
    } else {
      obs_all <- rbind(obs_all,data.frame(t(Y[ki,])))
    }
    
    yhat <- predict(model, Xtest) %>%
      data.frame() %>%
      setNames(colnames(Y)) 
    
    yhat_all <- rbind(yhat_all,yhat)

    # siteno and decade for joining
    site_no <- select(data[ki,],site_no,decade)
    site_no_all <- rbind(site_no_all,site_no)
    
  }
  
  average_mse <- data.frame(epoch = seq(1:epochs)) %>%
    mutate(val_mse = apply(mse_all, 2, mean))
  
    est_obs <- yhat_all %>%
      mutate(site_no=site_no_all$site_no,
             decade=site_no_all$decade) %>%
      gather(variable,nnet,-site_no,-decade) %>%
      mutate(obs = gather(obs_all,variable,obs)$obs)
  
  if(length(Ytest)==1){
    result <- list(average_mse=average_mse,
                   yhat_all=yhat_all,
                   obs_all=obs_all,
                   est_obs=est_obs,
                   time=proc.time() - ptm)
  }else{
    result <- list(average_mse=average_mse,
                   yhat_all=yhat_all,
                   obs_all=obs_all,
                   est_obs=est_obs,
                   #violations=violation_all,
                   time=proc.time() - ptm)
  }

  return(result)
  close(pb)
  
}




# area <- log10(data$tot_basin_area)
# null <- data.frame(L1=mean(Y[-ki,1]/area[-ki])*area[ki],
#                    T2=mean(Y[-ki,2]/area[-ki])*area[ki],
#                    T3=mean(Y[-ki,3]/area[-ki])*area[ki],
#                    T4=mean(Y[-ki,4]/area[-ki])*area[ki]) 
# 
# null_all <- rbind(null_all,null)

# custom loss function
# library(keras)
# y_true <- k_variable(value = array(c(3,6,8,9)))
# y_pred <- k_variable(value = array(c(4,5,3,7)))

# custom_loss <- function(y_true, y_pred){
#   var <- k_variable(value = array(y_pred))
#   n <- k_get_variable_shape(var)[[1]]
#   diff <- k_variable(value=array(y_pred[2:n]-y_pred[1:n-1]))
#   val <- k_get_value(k_sum(k_cast(k_less(diff,0),'float32')))
#   mse <- k_mean(k_square(y_true-y_pred)) 
#   error <- k_get_value(mse + val)
#   return(error)
# }
# 
# custom_loss <- function(y_true, y_pred){
#   n <- k_get_variable_shape(y_pred)[[1]]
#   diff <- k_eval(y_pred)[2:n]-k_eval(y_pred)[1:n-1]
#   val <- k_get_value(k_sum(k_cast(k_less(diff,0),'float32')))
#   mse <- k_mean(k_square(y_true-y_pred)) 
#   error <- k_get_value(mse + val)
#   return(error)
# }
