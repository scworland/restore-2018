
sw_dropout_uncertainty <- function(gage_all, iter=100) {
  
  # sourc k fold cv script
  
  sw_k_foldcv <- function(build_model,k=10,epochs=75,batch_size=25,X=X,Y=Y,data=d,loss_weights){
    
    library(keras)
    library(dplyr)
    library(tidyr)
    library(purrr)
    
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
      
      # early stopping callback
      # early_stopping <- callback_early_stopping(monitor = 'val_loss', patience = 2)
      
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
              verbose=0,
              callbacks=violation_count)
        
        violation_all <- rbind(violation_all,violation_count$val)
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
    
    mse_all <- data.frame(t(mse_all)) %>%
      set_names(paste0("kfold",1:k)) %>%
      mutate(epoch = seq(1:epochs)) %>%
      select(epoch,everything())
    
    est_obs <- yhat_all %>%
      mutate(site_no=site_no_all$site_no,
             decade=site_no_all$decade) %>%
      gather(variable,nnet,-site_no,-decade) %>%
      mutate(obs = gather(obs_all,variable,obs)$obs)
    
    if(length(Ytest)==1){
      result <- list(average_mse=average_mse,
                     mse_all=mse_all,
                     yhat_all=yhat_all,
                     obs_all=obs_all,
                     est_obs=est_obs,
                     time=proc.time() - ptm)
    }else{
      result <- list(average_mse=average_mse,
                     mse_all=mse_all,
                     yhat_all=yhat_all,
                     obs_all=obs_all,
                     est_obs=est_obs,
                     violations=violation_all,
                     time=proc.time() - ptm)
    }
    
    return(result)
    close(pb)
    
  }
  
  # start of dropout runs ----
  
  d <- gage_all
  
  f27 <- c("f0.02","f0.05","f0.1","f0.2","f0.5","f01","f02","f05",
           "f10","f20","f25","f30","f40","f50","f60","f70","f75",
           "f80","f90","f95","f98","f99","f99.5","f99.8","f99.9",
           "f99.95","f99.98")
  
  # grab basin area for later
  area <- d$basin_area
  
  # quantiles
  Y <- select(d,f27) %>%
    mutate_all(funs(.+0.001)) %>%
    mutate_all(funs(./area)) %>%
    mutate_all(funs(log10)) 
  
  X <- select(d,major:flood_storage, -basin_area) %>%
    mutate_all(funs(as.numeric(as.factor(.)))) %>%
    mutate_all(funs(as.vector(scale(.)))) %>%
    select_if(~!any(is.na(.)))  %>% 
    as.matrix()
  
  loss_weights <- c(rep(1,10),rep(1e-2,17))
  
  # Build model function for multiple outputs
  build_model <- function(){
    input <- layer_input(shape=dim(X)[2],name="basinchars")
    
    base_model <- input  %>%
      layer_dense(units = 28,activation="relu") %>%
      layer_dropout(rate=0.27) %>%
      layer_dense(units = 40,activation="relu") %>%
      layer_dropout(rate=0.70)  
    
    for(i in 1:dim(Y)[2]){
      y <- colnames(Y)[i]
      outstring <- paste0(
        sprintf("%s <- base_model %%>%%", y), 
        sprintf(" layer_dense(units = 1, activation='linear', name='%s')",y)
      )
      eval(parse(text=outstring))
    }
    
    Ylist <- paste0("list(",paste(colnames(Y),sep="",collapse=","),")")
    model <- keras_model(input,eval(parse(text=Ylist))) %>%
      compile(optimizer_rmsprop(lr = 0.0004),
              loss_weights = loss_weights,
              loss=keras::loss_mean_squared_error,
              metrics="mae")
    
    return(model)
  }
  
  cl <- makeCluster(5) # change
  registerDoParallel(cl)
  run <- foreach(i=1:iter) %dopar% {
    cv_results <- sw_k_foldcv(build_model,k=2,epochs=400,batch_size=100,Y=Y,X=X,data=d,loss_weights)
    cv_results$est_obs
  }
  stopCluster(cl)
  
  runs_combined <- run %>%
    purrr::reduce(left_join, by = c("site_no","decade","variable","obs")) %>%
    select(site_no,decade,variable,obs,everything()) %>%
    set_names(c("site_no","decade","variable","obs",paste0("run_",1:iter))) %>%
    gather(run,value,-site_no,-decade,-variable,-obs) %>%
    left_join(select(d,site_no,area=basin_area), by="site_no") %>%
    mutate(value = 10^(value)*area,
           obs = 10^(obs)*area,
            variable = as.numeric(substring(variable, 2))) #%>%
    # group_by(site_no,decade,variable) %>%
    # summarize(obs = mean(obs),
    #           mu = mean(value),
    #           min = min(value),
    #           max = max(value))
  
}