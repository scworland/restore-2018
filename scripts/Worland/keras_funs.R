
k_foldcv <- function(k=10,N,epochs=30,batch_size=250,data=d){
  
  #start the clock and progress bar
  ptm <- proc.time()
  pb <- txtProgressBar(min = 0, max = k, style = 3)
  
  k <- k
  indices <- sample(1:N)
  folds <- cut(indices, breaks = k, labels = FALSE)
  
  mse_all <- NULL
  obs_all <- NULL
  yhat_all <- NULL
  null_all <- NULL
  site_no_all <- NULL
  for (i in 1:k) {
    
    if(i %% 2==0) {
      # print progress every 2nd iteration
      setTxtProgressBar(pb, i)
    }
    
    ki <- which(folds == i)
    n <- length(ki)
    
    # training set
    Xtrain <- X[-ki,]
    Ytrain <- keras_array(list(Y[-ki,1],Y[-ki,2],Y[-ki,3],Y[-ki,4]))
    
    # validation set
    Xtest <- array_reshape(X[ki,],c(n,ncol(X)))
    Ytest <- keras_array(list(Y[ki,1],Y[ki,2],Y[ki,3],Y[ki,4]))
    
    # build and fit model
    model <- build_model()
    model_fit <- model %>% 
      fit(x=Xtrain,
          y=Ytrain, 
          epochs=epochs, 
          batch_size = batch_size,
          validation_data = list(Xtest, Ytest), 
          verbose=0)
    
    mse_all <- rbind(mse_all, model_fit$metrics$val_loss)
    
    average_mse <- data.frame(epoch = seq(1:epochs)) %>%
      mutate(val_mse = apply(mse_all, 2, mean))
    
    if(k < N){
      obs_all <- rbind(obs_all,data.frame(Y[ki,]))
    } else {
      obs_all <- rbind(obs_all,data.frame(t(Y[ki,])))
    }
    
    yhat <- predict(model, Xtest) %>%
      data.frame() %>%
      setNames(c("L1","T2","T3","T4")) 
    
    yhat_all <- rbind(yhat_all,yhat)
    
    area <- log10(data$tot_basin_area)
    null <- data.frame(L1=mean(Y[-ki,1]/area[-ki])*area[ki],
                       T2=mean(Y[-ki,2]/area[-ki])*area[ki],
                       T3=mean(Y[-ki,3]/area[-ki])*area[ki],
                       T4=mean(Y[-ki,4]/area[-ki])*area[ki]) 
    
    null_all <- rbind(null_all,null)

    # siteno and decade for joining
    site_no <- select(data[ki,],site_no,decade)
    site_no_all <- rbind(site_no_all,site_no)
  }
  
  est_obs <- yhat_all %>%
    mutate(site_no=site_no_all$site_no,
           decade=site_no_all$decade) %>%
    gather(lmom,neural_network,-site_no,-decade) %>%
    mutate(obs = gather(obs_all,lmom,obs)$obs,
           drainage_area = gather(null_all,lmom,null)$null) 
  
  result <- list(average_mse=average_mse,
                 yhat_all=yhat_all,
                 obs_all=obs_all,
                 null_all=null_all,
                 est_obs=est_obs,
                 time=proc.time() - ptm)
  
  return(result)
  close(pb)
  
}
