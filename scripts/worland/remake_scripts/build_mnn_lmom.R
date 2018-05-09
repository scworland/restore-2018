
sw_build_mnn_lmom <- function(gage_all) {
  
  d <- gage_all
  
  Y <- select(d,L1:T4) %>% 
    mutate(T2 = L2/L1,
           L1 = log10(L1)) %>%
    select(L1,T2,T3,T4) 
  
  X <- select(d,lon,lat,acc_hdens:statsgo) %>%
    mutate_at(vars(bedperm:statsgo),funs(as.numeric(as.factor(.)))) %>%
    mutate_all(funs(as.vector(scale(.)))) %>%
    select_if(~!any(is.na(.))) %>% 
    as.matrix()
  
  # Build model function for multiple outputs
  build_model <- function(){
    input <- layer_input(shape=dim(X)[2],name="basinchars")
    
    base_model <- input  %>%
      layer_dense(units = 40, activation='relu') 
    
    l1_pred <- base_model %>% 
      layer_dense(units = 1, activation='relu', name="l1") 
    
    t2_pred <- base_model %>% 
      layer_dense(units = 1, activation='relu', name="t2") 
    
    t3_pred <- base_model %>% 
      layer_dense(units = 1, activation='relu', name="t3") 
    
    t4_pred <- base_model %>% 
      layer_dense(units = 1, activation='relu', name="t4") 
    
    # pplo_pred <- base_model %>% 
    #   layer_dense(units = 1, activation='softplus', name="pplo") 
    
    model <- keras_model(input,list(l1_pred,t2_pred,t3_pred,t4_pred)) %>%
      compile(optimizer = "rmsprop",
              loss="mse",
              metrics="mae",
              loss_weights=c(0.2,0.5,0.5,1))
    
    return(model)
  }
  
  cv_results <- sw_k_foldcv(build_model,k=5,epochs=100,batch_size=50,Y=Y,X=X,data=d)
  
  # tail(cv_results$average_mse)
  # 
  # ggplot(cv_results$average_mse) +
  #   geom_line(aes(x = epoch, y = val_mse))
  # 
  # plots estimated vs observed
  est_obs <- cv_results$est_obs %>%
    # mutate(nnet = ifelse(variable=="L1",10^nnet,nnet),
    #        obs = ifelse(variable=="L1",10^obs,obs)) %>%
    gather(model,value,-variable,-obs,-site_no,-decade) %>%
    mutate(col=ifelse(obs>value,1,0))
  
  ggplot(est_obs) +
    geom_point(aes(obs,value,color=decade),alpha=0.7,shape=20) +
    scale_color_viridis_d() +
    geom_abline(slope=1,intercept=0,linetype="dashed",color="black") +
    facet_wrap(~variable,scales = "free",ncol=2) +
    labs(x="observed L-moment",y="estimated L-moment") +
    #ggtitle("10 fold CV") +
    theme_bw() +
    theme(legend.position = 'none')
  
  # subtitle="where drainage area = mean(log10(lmom[-k])/log10(area[-k]))*log10(area[k])"
  
  lmom_est <- cv_results$est_obs %>%
    mutate(nnet = ifelse(variable=="L1",10^nnet,nnet),
           obs = ifelse(variable=="L1",10^obs,obs)) %>%
    select(-obs) %>%
    spread(variable,nnet) %>%
    select(site_no,decade,L1,T2,T3,T4)
  
  # calculate quantiles from lmoms
  FF <- c(0.0002, 0.0005, 0.0010, 0.0020, 0.0050, 0.0100, 0.0200,
          0.0500, 0.1000, 0.2000, 0.2500, 0.3000, 0.4000, 0.5000, 
          0.6000, 0.7000, 0.7500, 0.8000, 0.9000, 0.9500, 0.9800,
          0.9900, 0.9950, 0.9980, 0.9990, 0.9995, 0.9998)
  
  lmom_list <- lmom_est %>%
    select(L1,T2,T3,T4) %>%
    as.matrix() %>%
    split(.,seq(nrow(.))) 
  
  # takes 10 minutes
  lmom_quants <- mclapply(lmom_list,sw_lmom2q,FF)
  
  # left pad and combine
  aep_all <- mclapply(lmom_quants,sw_left_na,num=27) %>%
    bind_rows() %>%
    t() %>%
    data.frame() %>%
    set_names(FF*100) %>%
    mutate(site_no=lmom_est$site_no,
           decade=lmom_est$decade) %>%
    gather(variable,aep,-site_no,-decade) %>%
    mutate(variable=as.numeric(variable),
           aep = ifelse(aep<0,0,aep)) 
  
  return(aep_all)
}



  