
sw_build_snn_quantile <- function(gage_all) {
  
  d <- gage_all
  
  f15 <- c("f0.02","f0.5","f05","f10","f20", "f30","f40","f50",
           "f60","f70","f80","f90","f95","f99.5","f99.98")
  
  Y <- select(d,f15) %>%
    mutate_all(funs(.+2)) %>%
    mutate_all(funs(log10))
  
  X <- select(d,major:flood_storage) %>%
    mutate_all(funs(as.numeric(as.factor(.)))) %>%
    mutate_all(funs(as.vector(scale(.)))) %>%
    select_if(~!any(is.na(.)))  %>% 
    as.matrix()
  
  # Build model function for multiple outputs
  build_model <- function(){
    input <- layer_input(shape=dim(X)[2],name="basinchars")
    
    base_model <- input  %>%
      layer_dense(units = 25,activation="relu") %>%
      layer_dropout(rate=0.558) %>%
      layer_dense(units = 49,activation="relu") %>%
      layer_dropout(rate=0.591) 
    
    output <- base_model %>%
      layer_dense(units = 1)
    
    model <- keras_model(input,output) %>%
      compile(optimizer_rmsprop(lr = 0.002536862),
              loss="mse",
              metrics="mae")
    
    return(model)
  }
  
  preds <- NULL
  for(j in 1:ncol(Y)){
    print(paste0("Processing quantile ",j, " out of ", ncol(Y)))
    y <- Y[,j]
    cv_results <- sw_k_foldcv(build_model,k=5,epochs=162,batch_size=294,Y=y,X=X,data=d)
    preds[[j]] <- cv_results$est_obs
  }
  
  preds_df <- bind_rows(preds) %>%
    mutate(obs = (10^obs)-2,
           nnet = round((10^nnet)-2,2),
           nnet = ifelse(nnet<0,0,nnet),
           variable = as.numeric(substring(variable, 2)))
  
  return(preds_df)
}


