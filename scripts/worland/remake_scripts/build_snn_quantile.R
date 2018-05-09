
sw_build_snn_quantile <- function(gage_all) {
  
  d <- gage_all
  
  Y <- select(d,f0.02:f99.98) %>%
    mutate_all(funs(.+2)) %>%
    mutate_all(funs(log10))
  
  X <- select(d,lon,lat,acc_hdens:statsgo) %>%
    mutate_at(vars(bedperm:statsgo),funs(as.numeric(as.factor(.)))) %>%
    mutate_all(funs(as.vector(scale(.)))) %>%
    select_if(~!any(is.na(.))) %>% 
    as.matrix()
  
  # Build model function for multiple outputs
  build_model <- function(){
    input <- layer_input(shape=dim(X)[2],name="basinchars")
    
    base_model <- input  %>%
      layer_dense(units = 40,activation="relu") %>%
      layer_dropout(rate = 0.1) %>%
      layer_dense(units = 30,activation="relu") %>%
      layer_dropout(rate = 0.1)
    
    output <- base_model %>%
      layer_dense(units = 1)
    
    model <- keras_model(input,output) %>%
      compile(optimizer = "rmsprop",
              loss="mse",
              metrics="mae")
    
    return(model)
  }
  
  preds <- NULL
  for(j in 1:ncol(Y)){
    print(paste0("Processing quantile ",j," out of 27"))
    y <- Y[,j]
    cv_results <- sw_k_foldcv(build_model,k=5,epochs=200,batch_size=30,Y=y,X=X,data=d)
    preds[[j]] <- cv_results$est_obs
  }
  
  preds_df <- bind_rows(preds) %>%
    mutate(obs = (10^obs)-2,
           nnet = round((10^nnet)-2,2),
           variable = as.numeric(substring(variable, 2)))
  
  return(preds_df)
}


