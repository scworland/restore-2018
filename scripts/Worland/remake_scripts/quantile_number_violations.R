
sw_count_violations <- function(gage_all){

  d <- gage_all
  
  # grab basin area for later
  area <- d$basin_area
  
  Y <- select(d,f0.02:f99.98) %>%
    mutate_all(funs(.+0.001)) %>%
    mutate_all(funs(./area)) %>%
    mutate_all(funs(log10))

  X <- select(d,major:flood_storage, -basin_area) %>%
    #select(-c(cor_vars$drop)) %>%
    mutate_all(funs(as.numeric(as.factor(.)))) %>%
    mutate_all(funs(as.vector(scale(.)))) %>%
    select_if(~!any(is.na(.)))  %>% 
    as.matrix()
  
  f27 <- c("f0.02","f0.05","f0.1","f0.2","f0.5","f01","f02","f05",
          "f10","f20","f25","f30","f40","f50","f60","f70","f75",
          "f80","f90","f95","f98","f99","f99.5","f99.8","f99.9",
          "f99.95","f99.98")
  
  f21 <- c("f0.02","f0.1","f0.5","f02","f05","f10","f20","f25",
          "f30","f40","f50","f60","f70","f75","f80","f90","f95",
          "f98","f99.5","f99.9","f99.98")
  
  f15 <- c("f0.02","f0.5","f05","f10","f20", "f30","f40","f50",
          "f60","f70","f80","f90","f95","f99.5","f99.98")
  
  
  f_all <- list(f27,f21,f15)
  
  viol_all <- NULL
  for (i in 1:length(f_all)) {
    
    print(paste0("Analyzing Scenario ",i," out of ", length(f_all)))
    
    # subset Y matrix
    Y2 <- Y[,f_all[[i]]]
    
    # Build model function for multiple outputs
    build_model <- function(){
      input <- layer_input(shape=dim(X)[2],name="basinchars")
      
      base_model <- input  %>%
        layer_dense(units = 28,activation="relu") %>%
        layer_dropout(rate=0.27) %>%
        layer_dense(units = 40,activation="relu") %>%
        layer_dropout(rate=0.70)  
      
      for(j in 1:dim(Y2)[2]){
        y <- colnames(Y2)[j]
        outstring <- paste0(
          sprintf("%s <- base_model %%>%%", y), 
          sprintf(" layer_dense(units = 1, activation='linear', name='%s')",y)
        )
        eval(parse(text=outstring))
      }
      
      Ylist <- paste0("list(",paste(colnames(Y2),sep="",collapse=","),")")
      model <- keras_model(input,eval(parse(text=Ylist))) %>%
        compile(optimizer_rmsprop(lr = 0.0004),
                loss="mse",
                metrics="mae")
      
      return(model)
    }
    
    # K-fold cross validation
    epochs <- 400
    cv_results <- sw_k_foldcv(build_model,k=10,epochs=epochs,batch_size=100,Y=Y2,X=X,data=d)
    
    viol_count <- data.frame(t(cv_results$violations)) %>%
      set_names(paste0("kfold",1:ncol(.))) %>%
      mutate(epoch=1:epochs) %>%
      gather(fold,value,-epoch) %>%
      mutate(quants_est = ncol(Y2))
    
    viol_all <- rbind(viol_all,viol_count)
  }
  
  return(viol_all)

  # viol_count <- viol_all %>%
  #   group_by(epoch,quants_est,quants_est) %>%
  #   summarize(min=min(value),
  #             max=max(value),
  #             mu=mean(value)) %>%
  #   mutate(quants_est = factor(quants_est,quants_est))
  # 
  # 
  # ggplot(viol_count) +
  #   geom_line(aes(x=epoch,y=mu,color=quants_est, group = quants_est)) +
  #   scale_color_viridis_d(begin=0.1,end=0.8,option="C") +
  #   coord_cartesian(xlim=c(0,150)) +
  #   theme_bw() +
  #   labs(color="#quantiles est", fill = "#quantiles est", y="number of violations")
}
