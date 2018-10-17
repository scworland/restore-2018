
sw_build_snn_quantile <- function(gage_all,cor_covars) {
  
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
  
  # remove 1 variable of highly correlated pairs
  # cor_vars <- distinct(cor_covars, drop)
  
  X <- select(d,major:flood_storage, -basin_area) %>%
    #select(-c(cor_vars$drop)) %>%
    mutate_all(funs(as.numeric(as.factor(.)))) %>%
    mutate_all(funs(as.vector(scale(.)))) %>%
    select_if(~!any(is.na(.)))  %>% 
    as.matrix()
  
  # Build model function for multiple outputs
  build_model <- function(){
    input <- layer_input(shape=dim(X)[2],name="basinchars")
    
    base_model <- input  %>%
      layer_dense(units = 20, activation="relu") %>%
      layer_dropout(rate=0.5) %>%
      layer_dense(units = 20, activation="relu") %>%
      layer_dropout(rate=0.5)  
    
    output <- base_model %>%
      layer_dense(units = 1)
    
    model <- keras_model(input,output) %>%
      compile(optimizer_rmsprop(lr = 0.0005),
              loss="mse",
              metrics="mae")
    
    return(model)
  }
  
  
  preds <- NULL
  for(j in 1:ncol(Y)){
    print(paste0("Processing quantile ",j," out of ", ncol(Y)))
    y <- Y[,j]
    
    cv_results <- sw_k_foldcv(build_model,k=2,epochs=400,batch_size=350,Y=y,X=X,data=d)
    
    pred_j <- cv_results$est_obs %>%
      mutate(obs = round(10^(obs)*area-0.001,2),
             nnet = round(10^(nnet)*area,2),
             #nnet = ifelse(nnet<0,0,nnet),
             variable = as.numeric(substring(variable, 2)))
    
    preds[[j]] <- pred_j
  }
  
  preds_df <- bind_rows(preds) 
  
  return(preds_df)
}


