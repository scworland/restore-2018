
sw_build_fdc_model <- function(gage_all, retrain=FALSE){
  
  if(isTRUE(retrain)){
    
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
    
    # train NN model
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
    
    # not very sensitive, just make it small
    loss_weights <- c(rep(1,10),rep(1e-2,17))
    
    Ylist <- paste0("list(",paste(colnames(Y),sep="",collapse=","),")")
    model <- keras_model(input,eval(parse(text=Ylist))) %>%
      compile(optimizer = optimizer_rmsprop(lr = 0.0004),
              loss_weights = loss_weights,
              loss="mse",
              metrics="mae")
    
    model_fit <- model %>%
      fit(x=X,
          y=Y,
          epochs=400,
          batch_size = 100,
          validation_split=0.1,
          verbose=0)
  }else{
    
    model <- load_model_hdf5("scripts/Worland/keras_training/models/fdc_model.hdf5")
    
  }
  
  return(model)
}