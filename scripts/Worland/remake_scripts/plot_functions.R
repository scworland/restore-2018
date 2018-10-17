
sw_plot_site_map <- function(simple_boundary,clean_sites,huc12s) {
  
  boundary <- simple_boundary
  gages <- clean_sites
  hucpp <- huc12s %>%
    select(lat,lon)
  
  regions = c("texas","alabama","florida","mississippi",
              "louisiana","georgia","arkansas","tennessee",
              "oklahoma","missouri","kentucky")
  
  states <- subset(map_data("state"),region %in% regions)
  
  cols <- c("gages"="dodgerblue", "pour_points"="grey50")
  
  ggplot() +
    geom_polygon(data=states, aes(x=long, y=lat, group=group),fill="grey85",
                 linetype="dotted",color="grey20",size=0.3) +
    geom_sf(data=boundary,fill="white",alpha=0.7,color="grey15") +
    coord_sf(crs = st_crs(boundary), datum = NA) +
    geom_point(data=hucpp, aes(lon,lat,color="pour_points"),shape=19,size=0.1) +
    geom_point(data=clean_sites, aes(lon,lat,color="gages"),shape=19, size=1.5) +
    scale_colour_manual(name="",values=cols) +
    #geom_point(data=hold, aes(dec_long_va,dec_lat_va), shape=21, size=2, fill="orange",color="red",alpha=0.6) +
    north(states, symbol = 3, scale = 0.10, location="topleft") +
    scalebar(states, dist = 250, dd2km = TRUE, model = 'WGS84', 
             st.size = 3, location="bottomleft") +
    theme_void() +
    theme(legend.position = c(0.72,0.25))
}


# plot period of record
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


# plot violations and error for each epoch
sw_plot_viol_error <- function(mnn_quantile_est,snn_quantile_est){
  
  snn_viol <- snn_quantile_est %>%
    group_by(site_no,decade) %>%
    summarize(viol = sum(cummax(nnet) != nnet))
  
  error <- mnn_quantile_est$mse_all %>%
    gather(kfold,"out-of-sample mean squared error",-epoch) 
  
  error_viol <- data.frame(t(mnn_quantile_est$violations)) %>%
    set_names(paste0("kfold",1:ncol(.))) %>%
    mutate(epoch=1:nrow(.)) %>%
    gather(kfold,"out-of-sample violation number",-epoch) %>%
    left_join(error, by=c("epoch","kfold")) %>%
    gather(type,value,-kfold,-epoch)
  
  ggplot(error_viol) +
    geom_line(aes(x=epoch,y=value,group=kfold),alpha=0.5) +
    # geom_hline(yintercept=sum(snn_viol$viol),linetype="dashed",color="grey40") +
    facet_wrap(~type, scales="free_y", ncol=1) +
    coord_cartesian(xlim=c(0,200)) +
    labs(y=NULL) +
    theme_bw()
}

# sw_plot_viol_error(mnn_quantile_est,snn_quantile_est)
  
# plot violation count for different number of quantiles
sw_plot_viol_count <- function(viol_count,snn_quantile_est) {
  
  snn_viol <- snn_quantile_est %>%
    group_by(site_no,decade) %>%
    summarize(viol = sum(cummax(nnet) != nnet))
  
  ggplot(viol_count) +
    # geom_hline(yintercept=sum(snn_viol$viol),linetype="dashed",color="grey40") +
    # geom_text(data = data.frame(),aes(x=200,y=sum(snn_viol$viol)+20,
    #                                   label='Number of violations from SNN model for 27 quantiles'),
    #           color="grey40") +
    geom_line(aes(x=epoch,y=value,color=as.character(quants_est),alpha=fold)) +
    scale_alpha_manual(values=rep(1,10),guide='none') +
    scale_color_viridis_d(begin=0.0,end=0.9,option="C") +
    coord_cartesian(xlim=c(0,350)) +
    theme_bw() +
    labs(color = "number of \nestimated \nquantiles", y="number of violations") +
    theme(legend.position = c(0.85, 0.65),
          legend.background = element_rect(fill="white",color="grey50"))
}

# plot example fdcs and overall error
sw_plot_fdcs_error <- function(mnn_quantile_est, snn_quantile_est, n=10){
  
  sites <- c('02320700-1960','07341500-1950','07351000-1960')  
  
  snn <- snn_quantile_est %>%
    select(site_no, decade, variable, snnet = nnet, -obs)
  
  fdcs <- mnn_quantile_est$est_obs %>% 
    left_join(snn,by=c("site_no","decade","variable")) %>%
    select(-area) %>%
    # group_by(site_no,decade) %>% 
    # nest() %>%
    # sample_n(n) %>%
    # #mutate(samp = map2(data, 10, sample_n)) %>% 
    # select(site_no, decade, data) %>%
    # unnest() %>%
    mutate(site_decade = paste0(site_no,"-",decade)) %>%
    filter(site_decade %in% sites) %>%
    select(-site_no,-decade) %>%
    rename(mnnet=nnet) %>%
    gather(type,value,-variable,-site_decade)
  
  p1 <- ggplot(fdcs) +
    geom_line(aes(x=variable,y=value,linetype=type,color=type)) +
    #geom_point(aes(x=variable,y=value,shape=type)) +
    scale_color_manual(values=c("blue","black","grey")) +
    scale_linetype_manual(values=c("dashed","solid","dashed")) +
    labs(y='Q',x='Non-exceedance Probability') +
    facet_wrap(~site_decade, scales="free_y",ncol=n) +
    scale_y_log10() +
    # scale_y_log10(labels = comma) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=1))
  
  snn <- snn_quantile_est %>%
    select(site_no, decade, variable, snnet = nnet, -obs) %>%
    group_by(site_no, decade) %>%
    mutate(snnet = sort(snnet)) %>%
    ungroup() %>%
    arrange(site_no,decade,variable)
  
  # plot average error
  error <- mnn_quantile_est$est_obs %>%
    left_join(snn,by=c("site_no","decade","variable")) %>%
    na.omit() %>%
    gather(model,value,-site_no,-decade,-variable,-obs,-area) %>%
    mutate(error=(abs(value-obs)/area)*1000) %>%
    group_by(variable, model) %>%
    summarize(low=quantile(error,0.05),
              med=quantile(error,0.5),
              high=quantile(error,0.95)) %>%
    ungroup() 
  
  p2 <- ggplot(error) +
    # geom_linerange(aes(variable,ymin=low,ymax=high,color=model),position=position_dodge(width=3)) +
    # geom_point(aes(variable,med,color=model),position=position_dodge(width=3)) +
    geom_ribbon(aes(variable,ymin=low,ymax=high,fill=model),alpha=0.3) +
    geom_line(aes(variable,low,group=model),linetype="solid",alpha=0.4) +
    geom_line(aes(variable,high,group=model),linetype="solid",alpha=0.4) +
    geom_line(aes(variable,med,color=model),linetype="dashed") +
    #geom_point(aes(variable,med,color=model),position=position_dodge(width=3)) +
    scale_fill_manual(values=c("blue","grey"), labels=c("mnnet","snnet")) +
    scale_color_manual(values=c("blue","black"), labels=c("mnnet","snnet")) +
    labs(x="Non-exceedance Probability",y="error (mm)") +
    scale_y_log10() +
    theme_classic()
  
  grid.arrange(p1,p2)
}

# sw_plot_fdcs_error(mnn_quantile_est, snn_quantile_est, n=3)

# plot maps of estimated values
sw_map_quantiles <- function(all_est_fdc, simple_boundary){
  
  regions = c("texas","alabama","florida","mississippi",
              "louisiana","georgia","arkansas","tennessee",
              "oklahoma","missouri","kentucky")
  
  states <- subset(map_data("state"),region %in% regions)
  
  map_points <- all_est_fdc %>%
    select(comid, lon, lat, decade, nep, q, area) %>%
    mutate(q = q+0.001) %>%
    filter(decade=="2000") 
  
  boundary <- simple_boundary

  ggplot() +
    geom_polygon(data=states, aes(x=long, y=lat, group=group),fill="grey85",
                 linetype="solid",color="grey20",size=0.3) +
    geom_sf(data=boundary,fill="white",alpha=0.7,color="grey15") +
    coord_sf(crs = st_crs(boundary), datum = NA) +
    geom_point(data=filter(map_points,nep==50.00), aes(lon,lat,color=log10(q)),alpha=0.7) +
    theme_void() +
    scale_color_viridis_c(option="C", name=expression(log10(m^{3})))

}

# sw_map_quantiles(all_est_fdc, simple_boundary)

# plot dropout uncertainty
sw_dropout_runs_plot <- function(dropout_uncertainty, n=3){
  
  #sites <- c('02294898−1990','02366000−1960','02377570−1990') 
  sites <- c('02320700-1960','07341500-1950','07351000-1960') 
  
  # rand_runs <- dropout_uncertainty %>%    
  #   group_by(site_no,decade) %>% 
  #   nest() %>%
  #   sample_n(n) %>%
  #   select(site_no, decade, data) %>%
  #   unnest() %>%
  #   mutate(site_decade = paste0(site_no,"-",decade)) %>%
  #   select(-site_no,-decade) 
  
  rand_runs <- dropout_uncertainty %>%    
    mutate(site_decade = paste0(site_no,"-",decade)) %>%
    filter(site_decade %in% sites) %>%
    select(-site_no,-decade) 
  
  p1 <- ggplot(rand_runs) +
    geom_line(aes(variable,value,group=run),alpha=0.05,color="blue") +
    geom_line(aes(variable,obs), size=0.5, linetype="dashed") +
    # geom_line(aes(variable,min),linetype="dashed") +
    # geom_line(aes(variable,max),linetype="dashed") +
    # geom_line(aes(variable,mu),linetype="dotted") +
    facet_wrap(~site_decade,scales="free_y") +
    labs(x="non-exceedance probability",y="Q") +
    theme_bw() +
    scale_y_log10()
  
  contains <- dropout_uncertainty %>%
    group_by(site_no,variable) %>%
    summarize(min=min(value),
              max=max(value),
              obs=mean(obs)) %>%
    ungroup() %>%
    mutate(contained = ifelse(obs>min & obs<max,1,0)) %>%
    group_by(variable) %>%
    add_count() %>%
    summarize(percent = sum(contained)/mean(n))
  
  p2 <- ggplot(contains) +
    geom_line(aes(variable,percent)) +
    geom_point(aes(variable,percent),shape=21,size=1.5) +
    labs(x="quantile",y="percent of observed contained") +
    coord_cartesian(ylim=c(0,1)) +
    theme_classic()
  
  grid.arrange(p1,p2)
}

# sw_dropout_runs_plot(dropout_uncertainty, n=3)

# plot lime results
sw_plot_lime <- function(lime_all, gage_all){
  
  explanation_all <- lime_all$explanation_all
  
  d <- gage_all
  
  f27 <- c("f0.02","f0.05","f0.1","f0.2","f0.5","f01","f02","f05",
           "f10","f20","f25","f30","f40","f50","f60","f70","f75",
           "f80","f90","f95","f98","f99","f99.5","f99.8","f99.9",
           "f99.95","f99.98")
  
  # grab basin area for later
  area <- d$basin_area
  
  quants_sd <- select(d,f27) %>%
    mutate_all(funs(.+0.001)) %>%
    mutate_all(funs(./area)) %>%
    mutate_all(funs(log10)) %>% 
    gather(quantile,value) %>%
    group_by(quantile) %>%
    summarize(stand_dev = sd(value))
  
  plot_lime <- explanation_all %>%
    select(-data) %>%
    add_count(feature) %>%
    left_join(quants_sd,by="quantile") %>%
    mutate(quantile = substring(quantile, 2),
           feature_weight = feature_weight/stand_dev) %>%
    #filter(n > 26000) %>%
    # filter(feature %in% c("bfi","ppt_mean","deciduous_forest")) %>%
    sample_n(20000) %>%
    group_by(feature) %>%
    mutate(muabs=abs(mean(feature_weight))) %>%
    ungroup() %>%
    mutate(feature=fct_reorder(feature, muabs,.desc = TRUE))
  
  ggplot(plot_lime) +
    #geom_line(aes(x=quantile,y=feature_weight, group=case),alpha=0.01) +
    geom_point(aes(x=quantile,y=feature_weight, group=case),
               position=position_dodge(width=0.2),alpha=0.01) +
    geom_hline(yintercept = 0, linetype="dashed") +
    facet_wrap(~feature, scales="free_y",ncol=1) +
    theme_bw() +
    labs(x="Non-exceedance Probability",
         y="feature weight/sd(obs quantile values)")
}

# sw_plot_lime(lime_all)

# plot bias in final layer
sw_plot_bias <- function(lime_all, gage_all){
  
  bias <- lime_all$bias
  
  d <- gage_all
  
  f27 <- c("f0.02","f0.05","f0.1","f0.2","f0.5","f01","f02","f05",
           "f10","f20","f25","f30","f40","f50","f60","f70","f75",
           "f80","f90","f95","f98","f99","f99.5","f99.8","f99.9",
           "f99.95","f99.98")
  
  # grab basin area for later
  area <- d$basin_area
  
  mu_bias <- select(d,f27) %>%
    mutate_all(funs(.+0.001)) %>%
    mutate_all(funs(./area)) %>%
    mutate_all(funs(log10)) %>% 
    gather(quantile,value) %>%
    group_by(quantile) %>%
    summarize(mu = mean(value)) %>%
    mutate(bias=bias,
           quantile = as.numeric(substring(quantile, 2))) %>%
    gather(variable,value,-quantile) 
  
  ggplot(mu_bias) + 
    geom_point(aes(quantile,value,shape=variable)) +
    scale_shape_manual(values=c(6,16)) +
    labs(x="Non-exceedance Probability") +
    theme_classic()
}

# map LIME results
sw_map_lime <- function(lime_all, gage_all, huc12_covariates){
  
  explanation_all <- lime_all$explanation_all
  
  regions = c("texas","alabama","florida","mississippi",
              "louisiana","georgia","arkansas","tennessee",
              "oklahoma","missouri","kentucky")
  
  states <- subset(map_data("state"),region %in% regions)
  
  map_points <- huc12_covariates %>%
    select(comid, lon=dec_long_va, lat=dec_lat_va, decade) %>%
    distinct(comid,decade,.keep_all=T) %>% # drop duplicate comids
    filter(decade==2000) %>%
    mutate(case = as.character(1:nrow(.)))
  
  d <- gage_all
  
  f27 <- c("f0.02","f0.05","f0.1","f0.2","f0.5","f01","f02","f05",
           "f10","f20","f25","f30","f40","f50","f60","f70","f75",
           "f80","f90","f95","f98","f99","f99.5","f99.8","f99.9",
           "f99.95","f99.98")
  
  # grab basin area for later
  area <- d$basin_area
  
  quants_sd <- select(d,f27) %>%
    mutate_all(funs(.+0.001)) %>%
    mutate_all(funs(./area)) %>%
    mutate_all(funs(log10)) %>% 
    gather(quantile,value) %>%
    group_by(quantile) %>%
    summarize(stand_dev = sd(value))
  
  map_lime <- explanation_all %>%
    select(case, quantile, feature, feature_weight) %>%
    left_join(quants_sd, by = "quantile") %>%
    mutate(weight = feature_weight/stand_dev) %>%
    left_join(map_points,by="case") %>%
    filter(quantile %in% c("f50"))
  
  p1 <- ggplot() +
    # geom_polygon(data=states, aes(x=long, y=lat, group=group),fill="grey85",
    #              linetype="solid",color="grey20",size=0.3) +
    #geom_sf(data=boundary,fill="white",alpha=0.7,color="grey15") +
    #coord_sf(crs = st_crs(boundary), datum = NA) +
    geom_point(data=filter(map_lime,feature=="hay_pasture"), 
               aes(lon,lat,color=weight),alpha=0.5) +
    theme_void() +
    ggtitle("Percent Hay and Pasture") +
    #labs(tag="A") +
    #scale_color_gradient(low="green4",high="gold") +
    scale_color_viridis_c(option="D")
  
  p2 <- ggplot() +
    # geom_polygon(data=states, aes(x=long, y=lat, group=group),fill="grey85",
    #              linetype="solid",color="grey20",size=0.3) +
    #geom_sf(data=boundary,fill="white",alpha=0.7,color="grey15") +
    #coord_sf(crs = st_crs(boundary), datum = NA) +
    geom_point(data=filter(map_lime,feature=="ppt_mean"), 
               aes(lon,lat,color=weight),alpha=0.5) +
    theme_void() +
    ggtitle("Mean PPT") +
    #labs(tag="B") +
    scale_color_gradient(low="magenta",high="yellow2")
  
  p3 <- ggplot() +
    # geom_polygon(data=states, aes(x=long, y=lat, group=group),fill="grey85",
    #              linetype="solid",color="grey20",size=0.3) +
    #geom_sf(data=boundary,fill="white",alpha=0.7,color="grey15") +
    # coord_sf(crs = st_crs(boundary), datum = NA) +
    geom_point(data=filter(map_lime,feature=="nid_storage"), 
               aes(lon,lat,color=weight),alpha=0.5) +
    theme_void() +
    ggtitle("NID storage") +
    #labs(tag="C") +
    scale_color_gradient(low="white",high="dodgerblue")
  
  p4 <- ggplot() +
    # geom_polygon(data=states, aes(x=long, y=lat, group=group),fill="grey85",
    #              linetype="solid",color="grey20",size=0.3) +
    #geom_sf(data=boundary,fill="white",alpha=0.7,color="grey15") +
    #coord_sf(crs = st_crs(boundary), datum = NA) +
    geom_point(data=filter(map_lime,feature=="temp_mean"), 
               aes(lon,lat,color=weight),alpha=0.5) +
    theme_void() +
    ggtitle("Mean temperature") +
    #labs(tag="D") +
    scale_color_gradient(low="slateblue3",high="yellow2")
  
  grid.arrange(p1,p2,p3,p4)
  
}

#sw_map_lime(lime_all, gage_all, huc12_covariates)
  
  
  