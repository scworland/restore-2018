

sw_build_mnn_quantile <- function(gage_all) {
  
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
  
  cv_results <- sw_k_foldcv(build_model,k=10,epochs=500,batch_size=100,Y=Y,X=X,data=d,loss_weights)
  
  cv_results$est_obs <- cv_results$est_obs %>%
    left_join(select(d,site_no,area=basin_area),by="site_no") %>%
    mutate(obs = round(10^(obs)*area-0.001,2),
           nnet = round(10^(nnet)*area,2),
           variable = as.numeric(substring(variable, 2))) %>%
    dplyr::distinct(site_no,decade,variable, .keep_all=TRUE)
  
  return(cv_results)
}

# viol_count <- data.frame(t(mnn_quantile_est$violations)) %>%
#   set_names(paste0("kfold",1:ncol(.))) %>%
#   mutate(epoch=1:nrow(.)) %>%
#   gather(kfold,value,-epoch) %>%
#   ggplot() +
#   geom_line(aes(epoch,value,color=kfold))

