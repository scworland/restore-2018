
sw_predict_fdcs <- function(gage_all,gage_covariates,huc12_covariates,mnn_quantile_est) {
  
  d <- gage_all
  
  Y <- select(d,f0.02:f99.98) %>%
    mutate_all(funs(.+2)) %>%
    mutate_all(funs(log10))
  
  X <- select(d,lon,lat,acc_hdens:statsgo,-perennial_ice_snow,-area_sqkm) %>%
    mutate_at(vars(bedperm:statsgo),funs(as.numeric(as.factor(.)))) %>%
    mutate_all(funs(as.vector(scale(.)))) %>%
    select_if(~!any(is.na(.))) %>% 
    as.matrix()
  
  # train NN model 
  input <- layer_input(shape=dim(X)[2],name="basinchars")
  
  base_model <- input  %>%
    layer_dense(units = 40,activation="relu") %>%
    layer_dropout(rate=0.1) %>%
    layer_dense(units = 30,activation="relu") %>%
    layer_dropout(rate=0.1) 
  
  for(i in 1:dim(Y)[2]){
    y <- colnames(Y)[i]
    outstring <- paste0(
      sprintf("%s <- base_model %%>%%", y), 
      sprintf(" layer_dense(units = 1, activation='relu', name='%s')",y)
    )
    eval(parse(text=outstring))
  }
  
  Ylist <- paste0("list(",paste(colnames(Y),sep="",collapse=","),")")
  model <- keras_model(input,eval(parse(text=Ylist))) %>%
    compile(optimizer = "rmsprop",
            loss="mse",
            metrics="mae")
  
  model_fit <- model %>% 
    fit(x=X,
        y=Y, 
        epochs=200, 
        batch_size = 30,
        validation_split=0.1, 
        verbose=0)
  
  # combine huc12 and all gage covariates
  all_covariates <- gage_covariates %>%
    bind_rows(huc12_covariates) %>% # bind rows
    select(comid,site_no,decade,everything()) %>%
    distinct(comid,decade,.keep_all=T) # drop duplicate comids
  
  # prepare covariates for model predictions
  all_x <- all_covariates %>% 
    select(lon,lat,acc_hdens:statsgo,-perennial_ice_snow,-area_sqkm,-decade) %>%
    mutate_at(vars(bedperm:statsgo),funs(as.numeric(as.factor(.)))) %>%
    mutate_all(funs(as.vector(scale(.)))) %>%
    #select_if(~!any(is.na(.))) %>% 
    as.matrix()
  
  cv_yhat <- mnn_quantile_est$est_obs %>%
    select(site_no,decade,f=variable,cv_q=nnet, obs_q=obs)
  
  # predict FDCs for every site
  predictions <- predict(model, all_x)
  
  yhat <- predictions %>%
    data.frame() %>%
    setNames(colnames(Y)) %>%
    mutate(comid=all_covariates$comid, # add comid
           site_no=all_covariates$site_no, # add site_no
           huc12=all_covariates$huc12, # add huc12s
           decade=as.character(all_covariates$decade), # add decade
           lon=all_covariates$lon, # add lon
           lat=all_covariates$lat) %>% # add lat
    gather(f,q,-comid,-site_no,-huc12,-decade,-lon,-lat) %>% # convert wide --> long
    mutate(f = as.numeric(substring(f, 2))) %>%
    left_join(cv_yhat, by = c("site_no","decade","f")) %>% # add cv estimates for gages
    mutate(q = ifelse(q<0,0,q), # remove negative q values
           q = (10^q)-2, # anti-log
           q = ifelse(!is.na(cv_q),cv_q,q), # replace with cv_q for gages
           q = round(q,2), # round predicted q
           obs_q = round(obs_q,2)) %>% # round observed
    select(-cv_q) %>% # remove cv q estimates
    group_by(comid,decade) %>%
    arrange(decade,comid,q) %>% # sort FDC for monotonicity
    #summarize(viol = sum(cummax(q)!=q)) # check violation number
    data.frame()
  
  return(yhat)
}

# write_feather(yhat,"data/gage/all_fdc_direct_est.feather")
