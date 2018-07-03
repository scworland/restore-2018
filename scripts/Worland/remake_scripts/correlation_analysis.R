
sw_plot_corr <- function(gage_all){
  
  # source("scripts/Worland/utils.R")
  
  # load data
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
  
  # correlation analysis with quantile
  corrr_analysis <- X %>%
    mutate(f0.02 = Y$f0.02,
           #f50 = Y$f50,
           f99.98 = Y$f99.98) %>%
    correlate() %>%
    focus(f0.02,f99.98) %>%
    mutate(mu = mean(c(f0.02,f99.98))) %>%
    rename(feature = rowname) %>%
    arrange(abs(mu)) %>%
    select(-mu) %>%
    mutate(feature = as.factor(feature)) %>%
    gather(quantile,value,-feature) %>%
    mutate(sign = ifelse(value>0,"pos","neg"))
  
  corrr_analysis %>%
    ggplot(aes(x = value, y = fct_reorder(feature, desc(value)))) +
    geom_segment(aes(xend = 0, yend = feature,color=sign)) +
    geom_vline(xintercept=0,linetype="dashed") +
    coord_cartesian(xlim=c(-0.8,0.8)) +
    geom_point(aes(shape=quantile,color=sign),fill="white") +
    scale_color_manual(values=c("blue","red"),guide=F) +
    scale_shape_manual(values=c(21,22)) +
    labs(x="correlation",y="features") +
    theme_bw()
  
  # # correlation EDA
  # x_cor <- X %>%
  #   #select(barren:woody_wetland) %>%
  #   correlate() %>% 
  #   focus(acc_basin_area) %>%
  #   arrange(desc(abs(acc_basin_area))) %>%
  #   na.omit()
  # 
  # rplot(x_cor)
}
