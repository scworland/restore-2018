

sw_build_mnn_quantile <- function(gage_all) {
  
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
      compile(optimizer_rmsprop(lr = 0.002536862),
              loss="mse",
              metrics="mae")
    
    return(model)
  }
  
  # K-fold cross validation
  # source("scripts/Worland/remake_scripts/keras_funs.R")
  
  cv_results <- sw_k_foldcv(build_model,k=10,epochs=162,batch_size=294,Y=Y,X=X,data=d)
  
  cv_results$est_obs <- cv_results$est_obs %>%
    mutate(obs = (10^obs)-2,
           nnet = round((10^nnet)-2,2),
           nnet = ifelse(nnet<0,0,nnet),
           variable = as.numeric(substring(variable, 2)))
  
  return(cv_results)
}

# tail(cv_results$average_mse)
# 
# ggplot(cv_results$average_mse) +
#   geom_line(aes(x = epoch, y = val_mse))
# 
# # fit metrics
# metrics <- cv_results$est_obs %>%
#   mutate(obs = (10^obs)-2,
#          nnet = round((10^nnet)-2,2),
#          variable = as.numeric(substring(variable, 2))) %>%
#   group_by(variable,decade) %>%
#   summarize(corr = cor(obs,nnet),
#             rmse = sqrt(mean((obs-nnet)^2))/mean(obs),
#             mpe = median(abs((nnet[obs>0]-obs[obs>0])/obs[obs>0]))*100) %>%
#   gather(metric,value,-variable,-decade)
#   
# ggplot(metrics) +
#   geom_line(aes(variable,value)) +
#   facet_grid(metric~decade,scales="free_y") +
#   labs(x = "exceedance probability") +
#   theme_bw()
# 
# # plots estimated vs observed
# est_obs <- cv_results$est_obs %>%
#   mutate(obs = (10^obs)-2,
#          nnet = round((10^nnet)-2,2),
#          variable = as.numeric(substring(variable, 2))) %>%
#   group_by(site_no,decade) %>%
#   mutate(nnet = sort(nnet)) %>%
#   # mutate(nnet = ifelse(nnet<0,0,nnet) %>%
#   #          ifelse(.==cummax(.),.,NA) %>%
#   #          na.approx(.,rule=2)) %>%
#          #loess_fit = predict(loess(nnet~variable,span=0.2))) %>%
#   gather(model,value,-site_no,-decade,-variable) %>%
#   group_by(site_no) %>%
#   sample_n_groups(6)
# 
# ggplot(est_obs) +
#   geom_line(aes(variable,value,linetype=model)) +
#   scale_linetype_manual(values=c("dashed","solid","dotted")) +
#   facet_wrap(site_no~decade,scales="free_y") +
#   labs(x="EP",y="Q") +
#   theme_bw() +
#   scale_y_log10()
# 
# # constrained GAM
# # hold <- cv_results$est_obs %>%
# #   rename(nn=model) %>%
# #   left_join(scales,by=c("variable","decade")) %>%
# #   mutate(nn = (nn*sigma)+mu-2) %>%
# #   select(-mu,-sigma) %>%
# #   filter(site_no == "02352500")  %>%
# #   mutate(variable = as.numeric(substring(variable, 2))) %>%
# #   group_by(decade,site_no) %>%
# #   do(augment(scam(nn ~ s(variable, k=14,bs="mpi"), data =.))) 
# 
# # Y <- select(d,decade,f0.02:f99.98) %>%
# #   mutate_all(funs(replace(., .==0, 0.001))) %>%
# #   mutate_all(funs(log10)) 
# 
# 
# hold <- data.frame(t(cv_results$violations)) %>%
#   set_names(paste0("kfold",1:ncol(.))) %>%
#   mutate(epoch=1:800) %>%
#   gather(fold,value,-epoch) 
# 
# ggplot(hold) +
#   geom_line(aes(epoch,value,group=fold),alpha=0.5) +
#   theme_bw() +
#   scale_y_log10()

