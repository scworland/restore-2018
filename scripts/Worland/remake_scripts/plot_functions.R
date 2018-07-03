
sw_plot_site_map <- function(simple_boundary,clean_sites,huc12s) {
  
  boundary <- simple_boundary
  gages <- clean_sites
  hucpp <- huc12s %>%
    select(lat,lon)
  
  regions = c("texas","alabama","florida","mississippi",
              "louisiana","georgia","arkansas","tennessee",
              "oklahoma","missouri","kentucky")
  
  states <- subset(map_data("state"),region %in% regions)
  
  ggplot() +
    geom_polygon(data=states, aes(x=long, y=lat, group=group),fill="grey85",
                 linetype="dotted",color="grey20",size=0.3) +
    geom_sf(data=boundary,fill="white",alpha=0.7,color="grey15") +
    coord_sf(crs = st_crs(boundary), datum = NA) +
    #geom_point(data=hucpp, aes(lon,lat),shape=1,size=0.1) +
    geom_point(data=clean_sites, aes(lon,lat),shape=4, size=1) +
    north(states, symbol = 3, scale = 0.10, location="topleft") +
    scalebar(states, dist = 250, dd2km = TRUE, model = 'WGS84', 
             st.size = 3, location="bottomleft") +
    theme_void() 
}

sw_plot_por <- function(dv_list) {
  
  yrs <- map_df(dv_list,`[`, c('site_no','Flow','wyear')) %>%
    group_by(site_no,wyear) %>%
    summarize(Q = log10(mean(Flow))) %>%
    ungroup() %>%
    spread(site_no,Q) %>%
    gather(site_no,Q,-wyear) %>%
    mutate(flow = ifelse(is.na(Q),"no","yes")) %>%
    group_by(site_no) %>%
    mutate(start=min(wyear[!is.na(Q)]),
           end=max(wyear[!is.na(Q)]),
           length=end-start) %>%
    arrange(start,length) %>%
    ungroup() %>%
    mutate(site_no = fct_reorder(site_no,start,.desc=T))
  
  ggplot(yrs) +
    geom_tile(aes(x=wyear, y=site_no, fill = flow),color="white", size=0.1) +
    # scale_fill_viridis_d(option="C",na.value="white") +
    scale_fill_manual(values=c("white","dodgerblue")) +
    labs(x="years",y="Individual gages", fill = "flow") +
    geom_vline(xintercept=1950,linetype="dashed",size=0.8) +
    theme_bw() +
    scale_x_continuous(expand=c(0,0)) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = c(.2, .25))
 
}

sw_plot_rmse <- function(mnn_quantile_est,gage_all){
  
  area <- gage_all %>%
    select(site_no,area=basin_area)
  
  rmse <- function(obs,est) {sqrt(mean((obs-est)^2,na.rm=T))}
  
  mnn_quantile_est$est_obs %>% 
    left_join(area, by='site_no') %>%
    mutate(obs = obs/area,
           nnet = nnet/area) %>%
    group_by(site_no,decade) %>% 
    summarize(r = rmse(obs,nnet)) %>%
    ungroup() %>%
    mutate(rbin = cut_number(r,n=4)) %>%
    ggplot() +
    geom_histogram(aes(r,fill=rbin),color='white',bins=200,alpha=0.8) + 
    scale_fill_viridis_d(begin=0.0,end=0.8,option="C") +
    labs(x="Drainage-area Normalized RMSE",fill='RMSE bins') +
    coord_cartesian(xlim=c(0,0.5)) +
    theme_bw() +
    theme(legend.position = c(0.85, 0.65),
          legend.background = element_rect(fill="white",color="grey50"))
  
}

sw_plot_fdcs <- function(mnn_quantile_est, gage_all){
  
  area <- gage_all %>%
    select(site_no,area=basin_area)
  
  rmse <- function(obs,est) {sqrt(mean((obs-est)^2,na.rm=T))}

    index <- mnn_quantile_est$est_obs %>% 
      left_join(area, by='site_no') %>%
      mutate(obs = obs/area,
             nnet = nnet/area) %>%
      group_by(site_no,decade) %>% 
      mutate(r = rmse(obs,nnet)) %>%
      ungroup() %>%
      mutate(rbin = cut_number(r,n=4)) %>%
      group_by(rbin) %>%
      nest() %>%
      mutate(samp = map2(data, 4, sample_n)) %>% 
      select(rbin, samp) %>%
      unnest() %>%
      mutate(site_decade = paste0(site_no,"-",decade))
  
    # one <- c('08175000-1970','02333500-1980','08068500-1960','08150800-2000')
    # two <- c('07266000-1960','07029500-1950','08166140-2000','02342933-1990')
    # three <- c('07281000-1970','08075000-1990','02421000-1960','02304500-1970')
    # four <- c('07287000-1950','08072730-1990','08114500-1950','02462500-2000')
    # 
    # sites <- c(one,two,three,four)
    
  fdc <- mnn_quantile_est$est_obs %>% 
    left_join(area, by='site_no') %>%
    mutate(obs = obs/area,
           nnet = nnet/area) %>%
    select(-area) %>%
    group_by(site_no,decade) %>% 
    mutate(r = rmse(obs,nnet)) %>%
    ungroup() %>%
    mutate(rbin = cut_number(r,n=4)) %>%
    mutate(site_decade = paste0(site_no,"-",decade)) %>%
    filter(site_decade %in% index$site_decade) %>%
    select(-site_no,-decade) %>%
    gather(type,value,-variable,-rbin,-site_decade,-r) %>%
    ungroup() %>%
    mutate(site_decade = fct_reorder(site_decade,r))
  
  ggplot(fdc) +
    geom_line(aes(x=variable,y=value,linetype=type,color=rbin)) +
    scale_color_viridis_d(begin=0.0,end=0.8,option="C") +
    scale_linetype_manual(values=c("dashed","solid")) +
    labs(color='RMSE',y='Q',x='Non-exceedance Probability',color='RMSE bins') +
    facet_wrap(~site_decade, ncol=4) +
    scale_y_log10(labels = comma) +
    theme_bw() +
    theme(legend.position = "top")
  

}

sw_plot_viol_count <- function(viol_count,snn_quantile_est) {
  
  snn_viol <- snn_quantile_est %>%
    group_by(site_no,decade) %>%
    summarize(viol = sum(cummax(nnet) != nnet))
  
  p <- viol_count %>%
    group_by(epoch,quants_est,quants_est) %>%
    summarize(min=min(value),
              max=max(value),
              mu=mean(value)) %>%
    mutate(quants_est = factor(quants_est,quants_est))
  
  
  ggplot(p) +
    geom_hline(yintercept=sum(snn_viol$viol),linetype="dashed",color="grey40") +
    geom_text(data = data.frame(),aes(x=45,y=6000,label='Number of violations from SNN model for 15 quantiles'),
              color="grey40") +
    geom_line(aes(x=epoch,y=mu,color=quants_est, group = quants_est)) +
    scale_color_viridis_d(begin=0.0,end=0.9,option="C") +
    coord_cartesian(xlim=c(0,100)) +
    theme_bw() +
    labs(color = "number of \nestimated \nquantiles", y="number of violations") +
    theme(legend.position = c(0.85, 0.65),
          legend.background = element_rect(fill="white",color="grey50"))
}

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
    select_if(~!any(is.na(.))) 
  
  # correlation analysis with quantile
  corrr_analysis <- X %>%
    mutate(f0.02 = Y$f0.02,
           #f50 = Y$f50,
           f99.98 = Y$f99.98) %>%
    correlate() %>%
    focus(f0.02,f99.98) %>%
    mutate(mu = mean(c(f0.02,f99.98))) %>%
    dplyr::rename(feature = rowname) %>%
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
    scale_color_viridis_d(begin=0.1,end=0.6,option="C",guide=F) +
    scale_shape_manual(values=c(21,22)) +
    labs(x="correlation",y="features") +
    theme_bw()
}

sw_plot_lime <- function(lime_all){
  
  plot_lime <- lime_all %>%
    select(-data) %>%
    add_count(feature) %>%
    mutate(quantile = substring(quantile, 2)) %>%
    #         as.numeric(.)/100) %>%
    filter(n > 5000) %>%
    group_by(feature,quantile) %>%
    summarize(min = min(feature_weight),
              mu = mean(feature_weight),
              max = max(feature_weight)) %>%
    ungroup() %>%
    group_by(feature) %>%
    mutate(muabs=abs(mean(mu))) %>%
    ungroup() %>%
    mutate(feature=fct_reorder(feature, muabs,.desc = TRUE))
  
  #length(unique(plot_lime$feature))
  
  ggplot(plot_lime) +
    geom_linerange(aes(x=quantile, ymin=min,ymax=max,color=feature)) +
    geom_point(aes(x=quantile, y=mu, color=feature),size=1) +
    scale_color_viridis_d(begin=0.0,end=0.8,option="C",direction=-1,guide=F) +
    facet_wrap(~feature, scales="free_y",ncol=2) + 
    geom_hline(yintercept=0,linetype="dashed") +
    labs(x='Non-exceedance Probability',y='feature weight') +
    scale_x_discrete(breaks=c("0.02","02","50","98","99.98"),
                     labels=c("01","02","50","98","99")) +
    #coord_cartesian(ylim=c(-0.1,0.65)) +
    theme_bw()
}