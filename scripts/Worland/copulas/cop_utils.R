
# function to return eastern huc4 data
east_data <- function(gage_all){
  
  result <- gage_all %>%
    mutate(huc4 = substr(huc12,1,4)) %>%
    filter(decade=="2000",
           nzero==0) %>%
    filter(huc4 %in% c("0316")) %>%
    st_as_sf(coords=c("lon","lat"),crs=4269)
  
  return(result)
  
  #loc = ifelse(huc4 == "0316","east","west")
}

# function to return western huc4 data
west_data <- function(gage_all){
  
  result <- gage_all %>%
    mutate(huc4 = substr(huc12,1,4)) %>%
    filter(decade=="2000",
           nzero==0) %>%
    filter(huc4 %in% c("1201","1202","1203","1204")) %>%
    st_as_sf(coords=c("lon","lat"),crs=4269) %>%
    filter(site_no != "08074020")
  
  return(result)
}

# function to return huc4 dvs
subset_dv <- function(gages,dv_list){
  
  result <- dv_list[unique(gages$site_no)] %>%
    bind_rows() %>%
    filter(year >= 2000 & year <= 2009) %>%
    select(site_no,date=Date,year,
           month,Q=Flow)
  
  if("08041500" %in% result$site_no){
    addition <- data.frame(site_no="08041500",date=as.Date("2007-09-30"),
                           year=2007,month=9,Q=137.5,stringsAsFactors = F)
    result <- bind_rows(result,addition) %>%
      arrange(site_no,date)
  }
  
  return(result)
}

# function to map study area
huc4_map <- function(restore_hucs,east_huc4,west_huc4){
  
  # col2 <- c("#c03728","#828585")
  col2 <- c("#c03728","#919c4c")
  
  # subset coordinates
  mobile_hucs <- restore_hucs %>%
    filter(substr(HUC_12,1,4) %in% "0316") %>%
    select(geometry) %>%
    st_union() %>%
    st_simplify(preserveTopology = TRUE, dTolerance = 0.01)
  
  nst_hucs <- restore_hucs %>%
    filter(substr(HUC_12,1,4) %in% c("1201","1202","1203","1204")) %>%
    select(geometry) %>%
    st_union() %>%
    st_simplify(preserveTopology = TRUE, dTolerance = 0.01)
  
  # subset rstates
  regions = c("texas","alabama","mississippi",
              "louisiana","arkansas","oklahoma")
  states <- subset(map_data("state"),region %in% regions)
  
  # make plot
  ggplot() +
    geom_polygon(data=states, aes(x=long, y=lat, group=group),fill="grey85",
                 linetype="solid",color="grey20",size=0.3) +
    geom_sf(data=mobile_hucs,fill="white",alpha=0.7,color="grey15") +
    geom_sf(data=nst_hucs,fill="white",alpha=0.7,color="grey15") +
    geom_sf(data=east_huc4, shape=4, color=col2[2],size=2,alpha=1) +
    geom_sf(data=west_huc4, shape=4, color=col2[1],size=2,alpha=1) +
    geom_text(aes(-87.3,35.5,label="Mobile-Tombigbee"),color=col2[2]) +
    geom_text(aes(-99,31,label="Galveston-Trinity"),color=col2[1]) +
    coord_sf(datum = NA) +
    north(states, symbol = 3, scale = 0.10, location="topleft") +
    scalebar(states, dist = 200, dd2km = TRUE, model = 'WGS84', 
             st.size = 3, location="bottomright") +
    theme_void() 
}

# function to plot time series
plot_ts <- function(east_dvs,west_dvs,gage_all,yr=NA){
  
  dvs <- bind_rows(east_dvs,west_dvs) %>%
    left_join(select(gage_all,site_no,huc12)) %>%
    mutate(huc4 = substr(huc12,1,4),
           loc = ifelse(huc4 == "0316","east","west")) %>%
    group_by(site_no) %>%
    mutate(nep = pobs(Q))
  
  if(is.na(yr)){
    ggplot(dvs) +
      geom_line(aes(x=date,y=nep,color=site_no),alpha=0.1) +
      scale_color_viridis_d(guide=FALSE) +
      facet_wrap(~loc, ncol=1) +
      labs(y="Non-exceedance probability") +
      theme_bw()
  }else{
    ggplot(filter(dvs, year==yr)) +
      geom_line(aes(x=date,y=nep,color=site_no),alpha=0.2) +
      scale_color_viridis_d(guide=FALSE) +
      facet_wrap(~loc, ncol=1) +
      labs(y="Non-exceedance probability") +
      theme_bw()
  }
}

# function to return matrix of Q
huc4_discharge <- function(huc4_dvs, logarithm=F){
  
  if(isTRUE(logarithm)){
    result <- huc4_dvs %>%
      select(date,site_no,Q) %>%
      mutate(Q = log10((Q/35.3147)+0.001)) %>%
      spread(site_no,Q) %>%
      select(-date) %>%
      as.matrix()
  }else{
    result <- huc4_dvs %>%
      select(date,site_no,Q) %>%
      mutate(Q=(Q/35.3147)+0.001) %>%
      spread(site_no,Q) %>%
      select(-date) %>%
      as.matrix()
  }
  
  return(result)
}

# function to calculate NSE
nse <- function(est,obs){
  return((1-(sum((est-obs)^2,na.rm=T)/sum((obs-mean(obs))^2,na.rm=T))))
} 

nse <- function(est,obs){
  return((1-(sum((log(est+0.001)-log(obs+0.001))^2,na.rm=T)/sum((log(obs+0.001)-mean(log(obs+0.001)))^2,na.rm=T))))
  }

# best case scenario "model"
model_best <- function(Q, R, U, Z, sigma, date, sites, coords, area){

  dmat <- distm(coords,fun=distGeo)
  
  # for loop for computation
  result_all <- NULL
  for(i in 1:ncol(Q)){
    
    print(paste0("Calculating for site ",i," out of ",ncol(Q)))
          
    # find the index of the nearest geographic neighbor
    nni <- sort(dmat[,i],index.return=TRUE)$ix[2]
    
    # find index of highest correlation
    ri <- rev(sort(sigma[,i],index.return=TRUE)$ix)[2]
    
    # direct prediction
    DAR_nearest_neighbor <- R[,nni]*area[i]
    DAR_highest_rho <- R[,ri]*area[i]
    QPPQ_nearest_neighbor <- quantile(R[,i],U[,nni])*area[i]
    QPPQ_highest_rho <- quantile(R[,i],U[,ri])*area[i]
    
    # pseudo bivariate copula prediction
    corr_weighted_runoff_nn <- (sigma[i,nni]*R[,nni])*area[i]
    corr_weighted_runoff_rho <- (sigma[i,ri]*R[,ri])*area[i]
    
    # bivariate copula prediction
    bivar_norm_cop_nn <- quantile(R[,i],pnorm(sigma[i,nni]*Z[,nni]))*area[i]
    bivar_norm_cop_rho <- quantile(R[,i],pnorm(sigma[i,ri]*Z[,ri]))*area[i]
    # var_est <- 1-(sigma[i,ri]^2)
    # bivariate <- quantile(target,pnorm(mean_est))
    # bivar_low <- quantile(target, pnorm(mean_est - (2*sqrt(var_est))))
    # bivar_high <- quantile(target, pnorm(mean_est + (2*(sqrt(var_est)))))
    
    # inverse distance weighted
    w <- 1/(dmat[-i,i]^2)
    IDW_runoff <- ((R[,-i]%*%w)/sum(w))*area[i]
    IDW_log_runoff <- 10^((log10(R[,-i])%*%w)/sum(w))*area[i]
    IDW_nep <- quantile(R[,i],(U[,-i]%*%w)/sum(w))*area[i]
    
    # pseudo multivariate copula prediction
    sigma_R <- cov(R, method="spearman")
    rx <- R[,-i]
    mux <- apply(rx,2,mean)
    sigma_xx_R <- sigma_R[-i,-i]
    sigma_yx_R <- sigma_R[i,-i]
    inv_sigma_xx_R <- solve(sigma_xx_R)
    ry <- mean(R[,i]) + t((sigma_yx_R%*%inv_sigma_xx_R)%*%(t(rx[,])-mux))
    corr_weighted_runoff <- ry*area[i]
    
    # multivariate copula prediction 
    zx <- Z[,-i]
    sigma_xx <- sigma[-i,-i]
    sigma_yx <- sigma[i,-i]
    inv_sigma_xx <- solve(sigma_xx)
    zy <- t((sigma_yx%*%inv_sigma_xx)%*%t(zx))
    multivar_norm_cop <- quantile(R[,i],pnorm(zy))*area[i]
    # var_y <- sigma[i,i]-(sigma_yx%*%inv_sigma_xx%*%t(sigma_yx)[,])
    # multivariate <- quantile(target,pnorm(zy))
    # multivar_low <- quantile(target, pnorm(zy - (2*sqrt(var_y)[[1]])))
    # multivar_high <- quantile(target, pnorm(zy + (2*sqrt(var_y)[[1]])))
    
    # mean of R, log(R), and U
    avg_runoff <- apply(R[,-i],1,mean)*area[i]
    avg_log_runoff <- 10^apply(log10(R[,-i]),1,mean)*area[i]
    avg_nep <- quantile(R[,i],apply(U[,-i],1,mean))*area[i]
    
    # combine estimates
    result <- data.frame(date = date,
                         target = sites[i],
                         donor_rho = sites[ri],
                         rho = sigma[i,ri],
                         observed = Q[,i],
                         DAR_nearest_neighbor,
                         DAR_highest_rho,
                         QPPQ_nearest_neighbor,
                         QPPQ_highest_rho,
                         corr_weighted_runoff_nn,
                         corr_weighted_runoff_rho,
                         bivar_norm_cop_nn,
                         bivar_norm_cop_rho,
                         IDW_runoff,
                         IDW_log_runoff,
                         IDW_nep,
                         corr_weighted_runoff,
                         multivar_norm_cop,
                         avg_runoff,
                         avg_log_runoff,
                         avg_nep,
                         Utarget = U[,i],
                         stringsAsFactors = F)
    
    result_all <- rbind(result_all,result)
    
  }
  
  rownames(result_all) <- NULL
  
  
  return(result_all)
}

# plot ridge plot of NSE
plot_best_nse <- function(estimates_best_e, estimates_best_w){
  
  by_method <- bind_rows(estimates_best_e, estimates_best_w) %>%
    mutate(location = c(rep("east",nrow(estimates_best_e)),rep("west",nrow(estimates_best_w)))) %>%
    select(-c(date, donor_rho, rho, Utarget)) %>%
    gather(model,value,-observed,-target,-location) %>%
    group_by(model,target,location) %>%
    summarize(NSE = nse(value,observed)) %>%
    ungroup() %>%
    mutate(location=ifelse(location=="east","Mobile-Tombigbee","Galveston-Trinity"),
           model = forcats::fct_reorder(model,NSE))
  
  col2 <- c("#c03728", "#919c4c")
  
  ggplot(by_method) +
    geom_density_ridges(aes(x=NSE,y=model,fill=location,color=location),
                        scale = 0.9, rel_min_height = 0.001,
                        jittered_points = TRUE, 
                        point_shape = 4,
                        position=position_points_jitter(width = 0.05, height = 0.1),
                        alpha=0.5) +
    facet_wrap(~location,ncol=2) +
    labs(y=NULL) +
    theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
    scale_fill_manual(values=col2,guide=FALSE) +
    scale_color_manual(values=col2,guide=FALSE) +
    coord_cartesian(xlim=c(-1,1)) +
    scale_x_continuous(expand = c(0.1, 0.1)) +
    scale_y_discrete(expand = c(0.1, 0.1)) 
  
}

# plot ts of best estimates
plot_ts_best_estimate <- function(estimates_best_e,estimates_best_w){
  
  # both sites have rho around 0.85
  
  # prepare eastern site
  single_site_e <- estimates_best_e %>%
    select(-c(donor_rho, rho, Utarget)) %>%
    filter(target=="02464360") %>%
    gather(model,value,-target,-date,-observed) %>%
    filter(date >= "2009-03-01" & date <="2009-09-01") %>%
    mutate(basin="Mobile-Tombigbee") 
  
  # prepare western site
  single_site_w <- estimates_best_w %>%
    select(-c(donor_rho, rho, Utarget)) %>%
    filter(target=="08066300") %>%
    gather(model,value,-target,-date,-observed) %>%
    filter(date >= "2009-03-01" & date <="2009-09-01") %>%
    mutate(basin = "Galveston-Trinity")
  
  combined_sites <- bind_rows(single_site_e, single_site_w) %>%
    filter(model %in% c("multivar_norm_cop","bivar_norm_cop_rho","DAR_nearest_neighbor","IDW_nep","QPPQ_highest_rho")) 
  
  nse_vals <- combined_sites %>%
    group_by(model,basin) %>%
    summarize(NSE=round(nse(value,observed),3))
  
  #cols <- c("#c03728", "#919c4c","#828585", "#efe1c6","#f5c04a")
  
  ggplot(combined_sites) +
    geom_line(aes(date,observed),color="grey30") +
    geom_line(aes(date,value,color=basin)) +
    scale_color_manual(values=c("#c03728", "#919c4c"),guide=FALSE) +
    geom_text(data=nse_vals,aes(x=as.Date("2009-08-01"),
                                y=850,label=paste0("NSE = ",NSE))) +
    facet_grid(model~basin) +
    scale_y_log10() +
    theme_bw(base_size=10)
  
}

# log-liklihoods for normal and t-dist
compare_log_lik <- function(Z, sigma){
  
  mu <- rep(0,ncol(Z))
  
  # log likelihood for mvn
  llmvn <- function(mu,sigma) {
    -sum(apply(Z,1,dmvn,mu,sigma,log=TRUE))
  }
  
  # log likelihood for mvt
  llmvt <- function(mu,sigma,df){
    -sum(apply(Z,1,dmvt,mu,sigma,df,log=TRUE))
  }
  
  # hack to find df
  dfs <- 1:30
  lls <- NULL
  for(i in 1:length(dfs)){
    lls[i] <- llmvt(mu, sigma, df=i)
  }
  
  # plot df vs ll
  p <- ggplot() +
    geom_point(aes(dfs,lls)) +
    geom_line(aes(dfs,lls)) +
    theme_bw()
  
  # find minimum df
  df_min <- dfs[which.min(lls)]
  
  # cacluate log-likelihods
  ll_t <- llmvt(mu, sigma, df=df_min)
  ll_n <- llmvn(mu, sigma)
  
  # return values
  result <- list(ll_t=ll_t,
                 ll_n=ll_n,
                 df=df_min,
                 df_plot=p)
  
  return(result)
  
}

# condtional variance for normal and t-dist
conditional_variance <- function(Q, R, U, Z, sigma, date, sites, area, df){
  result_all <- NULL
  for(i in 1:ncol(Q)){
    
    # multivariate copula prediction 
    zx <- Z[,-i]
    sigma_xx <- sigma[-i,-i]
    sigma_yx <- sigma[i,-i]
    inv_sigma_xx <- solve(sigma_xx)
    zy <- t((sigma_yx%*%inv_sigma_xx)%*%t(zx))
    mu <- unname(quantile(R[,i],pnorm(zy))*area[i])
    var_n <- sigma[i,i]-(sigma_yx%*%inv_sigma_xx%*%t(sigma_yx)[,])
    var_t <- diag((df + zx%*%inv_sigma_xx%*%t(zx))/(df + (ncol(Q)-1))*as.numeric(var_n))
    varn_low <- quantile(R[,i], pnorm(zy - (2*sqrt(var_n)[[1]])))*area[i]
    varn_high <- quantile(R[,i], pnorm(zy + (2*sqrt(var_n)[[1]])))*area[i]
    vart_low <- quantile(R[,i], pnorm(zy - (2*sqrt(var_t)[[1]])))*area[i]
    vart_high <- quantile(R[,i], pnorm(zy + (2*sqrt(var_t)[[1]])))*area[i]
    
    # combine estimates
    result <- data.frame(date = date,
                         target = sites[i],
                         observed = Q[,i],
                         mu,
                         var_n,
                         var_t,
                         varn_low,
                         varn_high,
                         vart_low,
                         vart_high,
                         Utarget = U[,i])
    
    result_all <- rbind(result_all,result)
    
    rownames(result_all) <- NULL
  }
  
  return(result_all)
}

# plot variances across Us
var_across_u <- function(var_compare){
  
  by_u <- var_compare %>%
    select(target, Utarget, norm_dist = var_n, t_dist = var_t) 
  
  cols <- c("norm-dist"="dodgerblue", "t-dist"="black")
    
  ggplot(by_u) +
    geom_line(aes(Utarget,t_dist,group=target,color="t-dist"),alpha=0.03) +
    geom_line(aes(Utarget,norm_dist,group=target,color="norm-dist")) +
    scale_colour_manual(name="distributions",values=cols) +
    labs(y="condtional z-score variance") +
    scale_y_log10() +
    theme_bw()
  
}

# plot ts with conditional variance
plot_ts_ci <- function(var_compare){
  
  low <- var_compare %>%
    filter(target=="02471001") %>%
    select(-c(var_n, var_t, target, Utarget, varn_high,vart_high)) %>%
    filter(date >= "2009-03-01" & date <="2009-09-01") %>%
    gather(dist,low,-mu,-observed,-date) %>%
    mutate(dist = ifelse(dist=="varn_low","norm-dist","t-dist"))
  
  high <- var_compare %>%
    filter(target=="02471001") %>%
    select(date, varn_high, vart_high) %>%
    filter(date >= "2009-03-01" & date <="2009-09-01") %>%
    gather(dist,high,-date) %>%
    mutate(dist = ifelse(dist=="varn_high","norm-dist","t-dist"))
  
  ci <- left_join(low, high, by=c("date","dist")) %>%
    mutate(out = ifelse(observed>high | observed<low ,observed,NA),
           dist = ifelse(dist=="norm-dist","normal-copula","t-copula"))
  
  ggplot(ci) +
    geom_ribbon(aes(x=date,ymin=low,ymax=high),fill="grey80") +
    geom_line(aes(date,observed), color="darkorchid1",linetype="solid",size=0.6) +
    facet_wrap(~dist,ncol=1) +
    geom_point(aes(date,out),color="black", fill="orange",size=1.5,shape=21) +
    #scale_y_log10() +
    labs(y=expression(m^3~s^{-1})) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
}

# check percentiles
check_percentiles <- function(var_compare_e,var_compare_w){
  
  percent_table_e <- var_compare_e %>%
    mutate(norm_dist = ifelse(observed>varn_high | observed<varn_low,0,1),
           t_dist = ifelse(observed>vart_high | observed<vart_low,0,1),
           U = cut_number(Utarget,20)) %>%
    group_by(U) %>%
    add_count() %>%
    summarize("norm-dist" = sum(norm_dist)/unique(n),
              "t-dist" = sum(t_dist/unique(n)),
              Umean = mean(Utarget)) %>%
    ungroup() %>%
    select("norm-dist","t-dist",Umean) %>%
    gather(type,value,-Umean) %>%
    mutate(location="Mobile-Tombigbee")
  
  percent_table_w <- var_compare_w %>%
    mutate(norm_dist = ifelse(observed>varn_high | observed<varn_low,0,1),
           t_dist = ifelse(observed>vart_high | observed<vart_low,0,1),
           U = cut_number(Utarget,20)) %>%
    group_by(U) %>%
    add_count() %>%
    summarize("norm-dist" = sum(norm_dist)/unique(n),
              "t-dist" = sum(t_dist/unique(n)),
              Umean = mean(Utarget)) %>%
    ungroup() %>%
    select("norm-dist","t-dist",Umean) %>%
    gather(type,value,-Umean) %>%
    mutate(location="Galveston-Trinity")
  
  percent_table <- bind_rows(percent_table_e,percent_table_w)
  
  cols <- c("#c03728", "#919c4c","#828585", "#efe1c6","#f5c04a")
  cols <- c("#828585", "#f5c04a")
  
  ggplot(percent_table,aes(Umean,value,color=type)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values=cols,name="",label=c("normal copula","t copua")) +
    geom_hline(yintercept=0.95, linetype="dashed") +
    labs(x="Non-exceedance probability",
         y=expression(paste("percent contained within ", mu%+-%2~sigma))) +
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(breaks=seq(0.0,1.0,0.2),labels=seq(0.0,1.0,0.2)) +
    facet_wrap(~location) +
    theme_bw() +
    theme(legend.position = c(0.9,0.2))
}

# compare NSE performance across rhos
compare_across_rhos <- function(Q, R, U, Z, sigma, date, sites, coords, area){
  
  dmat <- distm(coords,fun=distGeo)
  
  # for loop for computation
  result_all <- NULL
  for(i in 1:ncol(Q)){
    
    print(paste0("calculating pairwise NSE for site ",i," of ",ncol(Q)))
    
    # select target
    target <- Q[,i]
    
    for(j in 2:ncol(Q)){
      
      # sorted correlation vector -target
      ri_m <- rev(sort(sigma[,i],index.return=TRUE)$ix)[j:ncol(Q)]
      ri_s <- rev(sort(sigma[,i],index.return=TRUE)$ix)[j]
      
      # find the index of the nearest geographic neighbor
      nni <- sort(dmat[,i],index.return=TRUE)$ix[j]
      
      # direct prediction
      QPPQ_highest_rho <- quantile(R[,i],U[,ri_s])*area[i]
      
      # bivariate copula prediction
      bivar_norm_cop_rho <- quantile(R[,i],pnorm(sigma[i,ri_s]*Z[,ri_s]))*area[i]
      
      # weighted correlation
      sigma_R <- sigma
      rx <- R[,ri_m]
      if(length(ri_m)>1){
        mux <- apply(rx,2,mean)
        sigma_xx_R <- sigma_R[ri_m,ri_m]
        sigma_yx_R <- sigma_R[i,ri_m]
        inv_sigma_xx_R <- solve(sigma_xx_R)
        ry <- mean(R[,i]) + t((sigma_yx_R%*%inv_sigma_xx_R)%*%(t(rx[,])-mux))
        corr_weighted_runoff <- ry*area[i]
      }else{
        mux <- mean(rx)
        sigma_xx_R <- sigma_R[ri_m,ri_m]
        sigma_yx_R <- sigma_R[i,ri_m]
        inv_sigma_xx_R <- solve(sigma_xx_R)
        ry <- mean(R[,i]) + t((sigma_yx_R%*%inv_sigma_xx_R)%*%(t(rx)-mux))
        corr_weighted_runoff <- ry*area[i]
      }

      
      # multivariate copula prediction
      zx <- Z[,ri_m]
      sigma_xx <- sigma[ri_m,ri_m]
      sigma_yx <- sigma[i,ri_m]
      inv_sigma_xx <- solve(sigma_xx)
      zy <- t((sigma_yx%*%inv_sigma_xx)%*%t(zx))
      multivar_norm_cop <- quantile(R[,i],pnorm(zy))*area[i]
      
      # inverse distance weighted
      w <- 1/(dmat[ri_m,i]^2)
      if(length(ri_m)>1){
        IDW_nep <- quantile(R[,i],(U[,ri_m]%*%w)/sum(w))*area[i]
        IDW_log_runoff <- 10^((log10(R[,ri_m])%*%w)/sum(w))*area[i]
      }else{
        IDW_nep <- quantile(R[,i],U[,ri_m]*w)*area[i]
        IDW_log_runoff <- 10^((log10(R[,ri_m])*w))*area[i]
      }
      
      
      # combine estimates
      result <- data.frame(target = sites[i],
                           IDW_log_runoff=nse(IDW_log_runoff,target),
                           QPPQ_highest_rho=nse(QPPQ_highest_rho,target),
                           bivar_norm_cop_rho=nse(bivar_norm_cop_rho,target),
                           multivar_norm_cop=nse(multivar_norm_cop,target),
                           corr_weighted_runoff=nse(corr_weighted_runoff,target),
                           IDW_nep=nse(IDW_nep,target),
                           rho_rank=j-1,
                           max_rho = sigma[i,ri_s])
      
      result_all <- rbind(result_all,result)
    }
  }
  
  return(result_all)
}

# plot NSE vs rho for each method
plot_rho_compare <- function(rho_compare_e, rho_compare_w, by=3){
  
  ranks_e <- unique(rho_compare_e$rho_rank)
  n_ranks_e <- length(ranks_e)
  
  rhos_e <- rho_compare_e %>%
    mutate(rho_rank=as.factor(rho_rank)) %>%
    group_by(rho_rank) %>%
    mutate(rho=as.character(round(mean(max_rho),2))) %>%
    ungroup() %>%
    gather(model,NSE,-target,-rho_rank, -max_rho, -rho) %>%
    filter(rho_rank %in% c(seq(1,n_ranks_e,by=by),n_ranks_e)) %>%
    mutate(basin="Mobile-Tombigbee")
  
  ranks_w <- unique(rho_compare_w$rho_rank)
  n_ranks_w <- length(ranks_w)
  
  rhos_w <- rho_compare_w %>%
    mutate(rho_rank=as.factor(rho_rank)) %>%
    group_by(rho_rank) %>%
    mutate(rho=as.character(round(mean(max_rho),2))) %>%
    ungroup() %>%
    gather(model,NSE,-target,-rho_rank, -max_rho, -rho) %>%
    filter(rho_rank %in% c(seq(1,n_ranks_w,by=by),n_ranks_w)) %>%
    mutate(basin="Galveston-Trinity")
  
  cols <- c("#c03728", "#919c4c","#828585", "#efe1c6","#f5c04a","dodgerblue4")
  
  ggplot(bind_rows(rhos_e,rhos_w)) +
    geom_boxplot(aes(x=rho,y=NSE,fill=model),width=0.7,outlier.size = 0.3) +
    labs(x="mean maximum Spearman's correlation coefficient",
         y="NSE") +
    scale_fill_manual(values = cols, name="") +
    facet_wrap(~basin,ncol=1,scales="free_x") +
    theme_bw() +
    coord_cartesian(ylim=c(-1,1)) +
    theme(legend.position = "top")
  
}

# plot NSE vs rho for each method
plot_rho_compare2 <- function(rho_compare_e,rho_compare_w){
  
  cols <- c("#c03728", "#919c4c","#828585", "#efe1c6","#f5c04a","dodgerblue4")
  #cols <- c("#c03728", "#919c4c","#f5c04a")
  
  # comibne rhos
  rhos_e <- rho_compare_e %>%
    select(-rho_rank,-target) %>%
    gather(model,NSE,-max_rho) %>%
    mutate(basin="Mobile-Tombigbee")
  
  rhos_w <- rho_compare_w %>%
    select(-rho_rank,-target) %>%
    gather(model,NSE,-max_rho) %>%
    mutate(basin="Galveston-Trinity")
    
  ggplot(bind_rows(rhos_e,rhos_w)) +
    annotate("rect", xmin = 0.85, xmax = 0.9, ymin = -1, ymax = 1,
             fill = "dark grey", alpha = .5) +
    geom_point(aes(x=max_rho,y=NSE,color=basin),size=0.4,alpha=0.25) +
    labs(x="Spearman's correlation coefficient",
         y="NSE") +
    geom_hline(yintercept = 0, linetype="solid") +
    #geom_vline(xintercept = 0.85, linetype="dashed") +
    scale_color_manual(values = c("#c03728", "#919c4c"),name=NULL) +
    facet_wrap(~model,ncol=2,strip.position = "right") +
    guides(colour = guide_legend(override.aes = list(size=1.5,alpha=1))) +
    theme_bw(base_size=14) +
    coord_cartesian(ylim=c(-1,1)) +
    theme(legend.position = "top")
}


# compute IDW values
idw <- function(U,coords){
  
  for(i in 1:ncol(U)){
    target <- Q[,i]
    d <- distHaversine(coords[i,],coords[-i,])
    w <- 1/d^2
    uy <- (U[,-i]%*%w)/sum(w)
  }
  
}

# fit multiple copulas
fit_cops <- function(U,sigma){
  
  list_all <- list()
  for(i in 1:ncol(U)){
    
    # select target
    target <- U[,i]
    
    # sorted correlation vector -target
    ri_s <- rev(sort(sigma[,i],index.return=TRUE)$ix)[2]
    
    print(paste0("fitting normal copula for ",i," of ",ncol(U), " site pairs"))
    
    # fit normal copula
    mpl_n <- fitCopula(normalCopula(dim=2, dispstr="un"), U[,c(i,ri_s)], method="mpl")
    norm_loglik <- summary(mpl_n)$fitC@loglik
    
    print(paste0("fitting t copula for ",i," of ",ncol(U), " site pairs"))
    
    # fit t copula
    mpl_t <- fitCopula(tCopula(dim=2, dispstr="un"), U[,c(i,ri_s)], method="mpl")
    t_loglik <- summary(mpl_t)$fitC@loglik
    
    print(paste0("fitting gumbel copula for ",i," of ",ncol(U), " site pairs"))
    
    # fit gumbel copula
    mpl_g <- tryCatch({
      fitCopula(gumbelCopula(dim=2), U[,c(i,ri_s)], method="mpl")
    },
    error=function(e){
      message("no loglik due to bad optimization")
      return(NA)})
    if(!is.na(mpl_g)){
      gumb_loglik <- summary(mpl_g)$fitC@loglik
    }else{
      gumb_loglik <- NA
    }
    
    print(paste0("fitting clayton copula for ",i," of ",ncol(U), " site pairs"))
    
    # fit clayton copula
    mpl_c <- tryCatch({
      fitCopula(claytonCopula(dim=2), U[,c(i,ri_s)], method="mpl",start=1)
    },
    error=function(e){
      message("no loglik due to bad optimization")
      return(NA)})
    if(!is.na(mpl_c)){
      clay_loglik <- summary(mpl_c)$fitC@loglik
    }else{
      clay_loglik <- NA
    }
    
    print(paste0("fitting frank copula for ",i," of ",ncol(U), " site pairs"))
    
    # fit frank copula
    mpl_f <- tryCatch({
      fitCopula(frankCopula(dim=2),U[,c(i,ri_s)], method="mpl")
    },
    error=function(e){
      message("no loglik due to bad optimization")
      return(NA)})
    if(!is.na(mpl_f)){
      frank_loglik <- summary(mpl_f)$fitC@loglik
    }else{
      frank_loglik <- NA
    }
    
    list_all[[i]] <- data.frame(site1 = colnames(U)[i],
                                site2 = colnames(U)[ri_s],
                                normal = norm_loglik,
                                t = t_loglik,
                                gumbel = gumb_loglik,
                                clayton = clay_loglik,
                                frank = frank_loglik)
  }
  
  df_all <- bind_rows(list_all)
  return(df_all)
}

# plot logliks
plot_logliks <- function(cops_e,cops_w){
  
  all_logliks <- bind_rows(cops_e,cops_w) %>%
    mutate(location=c(rep("Mobile-Tombigbee",nrow(cops_e)),
                      rep("Galveston-Trinity",nrow(cops_w)))) %>%
    gather(copula,loglik,-location,-site1,-site2) %>%
    group_by(site1,site2) %>%
    mutate(loglik=loglik/max(loglik))
  
  cols <- c("#c03728", "#919c4c","#828585", "#efe1c6","#f5c04a")

  ggplot(all_logliks) +
    geom_density_ridges(aes(x=loglik,y=location,fill=copula),
                        scale = 0.9, rel_min_height = 0.001,alpha=0.5,
                        size=0.1) +
    labs(y=NULL,x="normalize log-likelihood") +
    theme_ridges(grid = FALSE, font_size = 9, 
                 center_axis_labels = TRUE) +
    #scale_color_manual(values=cols) +
    scale_fill_manual(values=cols) +
    guides(fill=guide_legend(ncol=2)) +
    guides(color=guide_legend(ncol=2)) +
    coord_cartesian(xlim=c(0.4,1.1)) +
    #scale_x_continuous(expand = c(0.1, 0.1)) +
    scale_y_discrete(expand = c(0.1, 0.1)) +
    theme(legend.position = c(0.05,0.40))
}

# plot variance of different copulas
plot_cop_vars <- function(cops_w,U_w,R_w,area){
  
  gumbel_area <- area[west_huc4$site_no==gumbel_sites$value]
  
  # subset sites
  gumbel_sites <- cops_w %>%
    mutate(diff = gumbel-normal) %>%
    filter(diff == max(diff)) %>%
    select(site1,site2) %>%
    slice(1) %>%
    gather()
  
  # 1828:2040 = "2005-01-01" to "2005-08-01"
  
  # subset runoff
  R_gumbel <- R_w[1828:2040,gumbel_sites$value]
  
  # subset psuedo observations
  U_gumbel <- U_w[1828:2040,gumbel_sites$value]
  u1 <- U_gumbel[,1]
  u2 <- U_gumbel[,2]
  
  # find parameters 
  gumbel_pars <- BiCopSelect(u1,u2,familyset=4)
  normal_pars <- BiCopSelect(u1,u2,familyset=1)
  
  # fit gumbel and normal copula
  fit_gumbel <-  BiCop(family=4, gumbel_pars$par, tau = NULL, check.pars = TRUE)
  fit_normal <-  BiCop(family=1, normal_pars$par, tau = NULL, check.pars = TRUE)
  
  #BiCopfunc1(0.3, 0.1, fit_gumbel)
  
  # conditional gumbel simulation
  gumbel_sims <- matrix(data=NA,nrow=length(u1),ncol=100)
  for(i in 1:length(u1)){
    gumbel_sims[i,] <- BiCopCondSim(100, u1[i], 1, fit_gumbel)
  }
  
  # conditional gumbel simulation
  normal_sims <- matrix(data=NA,nrow=length(u1),ncol=100)
  for(i in 1:length(u1)){
    normal_sims[i,] <- BiCopCondSim(100, u1[i], 1, fit_normal)
  }
  
  gumbel_sims <- data.frame(gumbel_sims) %>%
    set_names(paste0("run",1:100)) %>%
    mutate(date = date[1828:2040],
           observed = u2,
           area = gumbel_area[2]) %>%
    gather(run,value,-observed,-date,-area) %>%
    mutate(Qest = quantile(R_gumbel[,2],value)*area,
           Qobs = R_gumbel[,2]*area,
           copula = "gumbel")
  
  normal_sims <- data.frame(normal_sims) %>%
    set_names(paste0("run",1:100)) %>%
    mutate(date = date[1828:2040],
           observed = u2,
           area = gumbel_area[2]) %>%
    gather(run,value,-observed,-date,-area) %>%
    mutate(Qest = quantile(R_gumbel[,2],value)*area,
           Qobs = R_gumbel[,2]*area,
           copula = "normal")
  
  result <- bind_rows(gumbel_sims,normal_sims)
  
  ggplot(result) + 
    geom_line(aes(date,Qest,group=run),size=0.1,alpha=0.2) +
    geom_line(aes(date,Qobs),color="dodgerblue") +
    facet_wrap(~copula,ncol=1) +
    scale_y_log10() +
    theme_bw()
  
  normal_mu <- normal_sigma <- NULL
  for(i in 1:nrow(U_w)){
    sim <- BiCopCondSim(1000, u1[i], 1, fit_normal)
    normal_mu[i] <- mean(sim)
    normal_sigma[i] <- sd(sim)
  }
  
  # conditional gumbel simulation
  gumbel_mu <- gumbel_sigma <- NULL
  for(i in 1:nrow(U_w)){
    sim <- BiCopCondSim(1000, u1[i], 1, fit_gumbel)
    gumbel_mu[i] <- mean(sim)
    gumbel_sigma[i] <- sd(sim)
  }
  
  normal_mu <- normal_sigma <- NULL
  for(i in 1:nrow(U_w)){
    sim <- BiCopCondSim(1000, u1[i], 1, fit_normal)
    normal_mu[i] <- mean(sim)
    normal_sigma[i] <- sd(sim)
  }
  
  gumbel_result <- data.frame(date = date,
                              u = u2, 
                              area = gumbel_area[2],
                              R = R_gumbel[,2],
                              mu = gumbel_mu, 
                              sigma = gumbel_sigma,
                              copula="gumbel")
  
  normal_result <- data.frame(date = date,
                              u = u2, 
                              area = gumbel_area[2],
                              R = R_gumbel[,2],
                              mu = normal_mu, 
                              sigma = normal_sigma,
                              copula="normal")
  
  result <- bind_rows(gumbel_result,normal_result) %>%
    filter(date>"2005-01-01" & date<"2005-08-01") %>%
    mutate(Q = R*area,
           av = quantile(R,mu)*area,
           low = ifelse(mu-(2*sigma)<0,0,mu-(2*sigma)),
           high = ifelse(mu+(2*sigma)>1,1,mu+(2*sigma)),
           low = quantile(R,low)*area,
           high = quantile(R,high)*area) %>%
    select(date,copula,Q,av,low,high) 
  
  ggplot(result) +
    geom_ribbon(aes(date,ymin=low,ymax=high),fill="grey") +
    geom_line(aes(date,Q),color="blue") +
    facet_wrap(~copula,ncol=1) +
    theme_bw()
  
  
}

# compute tail dependence
tail_dep <- function(U_e,U_w,p=0.01,pll=pll){
  
  meth="Schmid.Schmidt"
  
  col2 <- c("#c03728", "#919c4c")
  
  # eastern tail dependence
  le <- tril(fitLambda(U_e, p = p,method=meth),-1)
  ue <- triu(fitLambda(U_e, p = p, method=meth, lower.tail = FALSE),1)
  lue <- as.data.frame(as.matrix(ue + le)) %>%
    set_names(colnames(U_e)) %>%
    mutate(site1=colnames(U_e)) %>%
    gather(site2,lambda,-site1) %>%
    filter(site1!=site2) %>%
    mutate(upper="upper",
           lower="lower")
  
  p_lue <- ggplot(data = lue, aes(x=site1, y=site2, fill=lambda)) + 
    facet_grid(lower~upper) +
    geom_tile(color="grey") +
    scale_fill_gradient(low="white",high=col2[2],
                        name=expression(lambda)) +
    ggtitle("Mobile-Tombigbee") +
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1, size=6),
          axis.text.y = element_text(size=6),
          axis.title = element_blank(),
          panel.background = element_rect(fill = 'grey65'))
  
  # western tail dependence
  lw <- tril(fitLambda(U_w, p = p,method=meth),-1)
  uw <- triu(fitLambda(U_w, p = p, method=meth, lower.tail = FALSE),1)
  luw <- as.data.frame(as.matrix(uw + lw)) %>%
    set_names(colnames(U_w)) %>%
    mutate(site1=colnames(U_w)) %>%
    gather(site2,lambda,-site1) %>%
    filter(site1!=site2) %>%
    mutate(upper="upper",
           lower="lower")
  
  p_luw <- ggplot(data = luw, aes(x=site1, y=site2, fill=lambda)) + 
    facet_grid(lower~upper) +
    geom_tile(color="grey") +
    scale_fill_gradient(low="white",high=col2[1],
                        name=expression(lambda), 
                        limits=c(0,0.85)) +
    ggtitle("Galveston-Trinity") +
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1, size=6),
          axis.text.y = element_text(size=6),
          axis.title = element_blank(),
          panel.background = element_rect(fill = 'grey65'))
  
  # select sites
  upper_dep_w <- c("08057000","08049500")
  lower_dep_w <- c("08068000","08062500")
  upper_dep_e <- c("02457595","02458600")
  lower_dep_e <- c("02443500","02469800")
  
  # for nonexceedance probs
  upper_Us_e <- data.frame(U_e[,upper_dep_e]) %>% 
    set_names(c("site1","site2")) %>%
    mutate(type="upper tail",
           location="Mobile-Tombigbee",
           sitea = upper_dep_e[1],
           siteb = upper_dep_e[2])
  
  lower_Us_e <- data.frame(U_e[,lower_dep_e]) %>% 
    set_names(c("site1","site2")) %>%
    mutate(type="lower tail",
           location="Mobile-Tombigbee",
           sitea = lower_dep_e[1],
           siteb = lower_dep_e[2])
  
  upper_Us_w <- data.frame(U_w[,upper_dep_w]) %>% 
    set_names(c("site1","site2")) %>%
    mutate(type="upper tail",
           location="Galveston-Trinity",
           sitea = upper_dep_w[1],
           siteb = upper_dep_w[2])
  
  lower_Us_w <- data.frame(U_w[,lower_dep_w]) %>% 
    set_names(c("site1","site2")) %>%
    mutate(type="lower tail",
           location="Galveston-Trinity",
           sitea = lower_dep_w[1],
           siteb = lower_dep_w[2])
  
  dep_all <- bind_rows(upper_Us_w,lower_Us_w,upper_Us_e,lower_Us_e) 
  
  site_names <- dep_all %>%
    distinct(sitea,siteb,location,type) %>%
    gather(site,name,-location,-type) %>%
    mutate(x=c(rep(0.8,4),rep(0.05,4)),
           y=c(rep(0.05,4),rep(0.8,4)))
  
  dep_Us <- ggplot(dep_all) +
    geom_point(aes(site1,site2,color=location),alpha=0.05,size=0.5) +
    scale_color_manual(values=col2,guide=FALSE) +
    facet_grid(type~location) +
    labs(x="",y="") +
    theme_bw()
  
  # # for flow
  upper_Qs_e <- data.frame(Q_e[,upper_dep_e]) %>%
    mutate(date=date) %>%
    filter(date >= "2000-04-01" & date <="2000-11-01") %>%
    gather(site,Q,-date) %>%
    mutate(site=str_sub(site,2)) %>%
    mutate(type="upper tail",
           location="Mobile-Tombigbee")

  lower_Qs_e <- data.frame(Q_e[,lower_dep_e]) %>%
    mutate(date=date) %>%
    filter(date >= "2000-04-01" & date <="2000-11-01") %>%
    gather(site,Q,-date) %>%
    mutate(site=str_sub(site,2)) %>%
    mutate(type="lower tail",
           location="Mobile-Tombigbee")
  
  upper_Qs_w <- data.frame(Q_w[,upper_dep_w]) %>%
    mutate(date=date) %>%
    filter(date >= "2000-04-01" & date <="2000-11-01") %>%
    gather(site,Q,-date) %>%
    mutate(site=str_sub(site,2)) %>%
    mutate(type="upper tail",
           location="Galveston-Trinity")
  
  lower_Qs_w <- data.frame(Q_w[,lower_dep_w]) %>%
    mutate(date=date) %>%
    filter(date >= "2000-04-01" & date <="2000-11-01") %>%
    gather(site,Q,-date) %>%
    mutate(site=str_sub(site,2)) %>%
    mutate(type="lower tail",
           location="Galveston-Trinity")
  
  dep_Q_all <- bind_rows(upper_Qs_w,lower_Qs_w,upper_Qs_e,lower_Qs_e) 

  dep_qs <- ggplot(dep_Q_all) +
    geom_line(aes(date,Q,color=location,group=site),alpha=0.7) +
    scale_y_log10() +
    scale_color_manual(values=col2,guide=FALSE) +
    facet_grid(type~location) +
    labs(x="",y="") +
    theme_bw()
  
  
  ggdraw() +
    draw_plot(p_luw, x = 0, y = 0.5, width = 0.5, height = 0.5) +
    draw_plot(p_lue, x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
    draw_plot(dep_Us, x = 0, y = 0, width = 0.41, height = 0.45) +
    draw_plot(dep_qs, x = 0.5, y = 0, width = 0.41, height = 0.45) +
    draw_plot_label(label = c("A.","B.","C."), size = 13,
                    x = c(0,0,0.5), y = c(1, 0.5, 0.5))
  
  
}

# plot simulated copulas
plot_sim_cops <- function(N=1000,alpha=0.1){
  
  # normal copula
  norm_cop <- normalCopula(param=0.8,dim=2) 
  nc_mv <- mvdc(norm_cop,margins=c("norm","norm"),
                paramMargins=list(list(0,2), list(0,5)))
  
  nc_x <- data.frame(rMvdc(N, nc_mv)) %>%
    set_names(c("x1","x2")) %>%
    mutate(copula="Normal")
  
  # t copula
  t_cop <- tCopula(param=0.8,dim=2,df = 1) 
  tc_mv <- mvdc(t_cop,margins=c("norm","norm"),
                paramMargins=list(list(0,2), list(0,5)))
  
  tc_x <- data.frame(rMvdc(N, tc_mv)) %>%
    set_names(c("x1","x2")) %>%
    mutate(copula="Student-t")
  
  # gumbel copula
  gumbel_cop <- gumbelCopula(param=5, dim=2) 
  gc_mv <- mvdc(gumbel_cop,margins=c("norm","norm"),
             paramMargins=list(list(0,2), list(0,5)))
  
  gc_x <- data.frame(rMvdc(N, gc_mv)) %>%
    set_names(c("x1","x2")) %>%
    mutate(copula="Gumbel")
  
  # clayton copula
  clayton_cop <- claytonCopula(param=5, dim=2) 
  cc_mv <- mvdc(clayton_cop,margins=c("norm","norm"),
                paramMargins=list(list(0,2), list(0,5)))
  
  cc_x <- data.frame(rMvdc(N, cc_mv)) %>%
    set_names(c("x1","x2")) %>%
    mutate(copula="Clayton") 
  
  # combine
  levels = c("Normal","Student-t","Gumbel","Clayton")
  df <- bind_rows(nc_x,tc_x,gc_x,cc_x) %>%
    mutate(copula = fct_relevel(copula,levels))
  
  ggplot(df) +
    geom_point(aes(x1,x2),alpha=alpha) +
    facet_wrap(~copula,ncol=4) +
    theme_bw(base_size=14) +
    labs(x=expression(x[1]), y=expression(x[2])) +
    theme(axis.text=element_blank(),
          axis.ticks=element_blank())
}


# estimated FDC and correlation matrix
model_fdc_est <- function(Q, R, U, Z, sigma, date, sites, coords, area, est_fdcs){
  
  dmat <- distm(coords,fun=distGeo)
  
  # for loop for computation
  result_all <- NULL
  for(i in 1:ncol(Q)){
    
    print(paste0("Calculating for site ",i," out of ",ncol(Q)))
    
    # subset estimated fdc
    fdc_est <- est_fdcs %>% 
      filter(site_no==sites[i] & decade=="2000") %>%
      mutate(nep=nep/100,
             cv_r=cv_q/area)
    
    # function for nep ---> Q
    nn_quantile <- approxfun(fdc_est$nep,fdc_est$cv_r)
    
    # find the index of the nearest geographic neighbor
    nni <- sort(dmat[,i],index.return=TRUE)$ix[2]
    
    # find index of highest correlation
    ri <- rev(sort(sigma[,i],index.return=TRUE)$ix)[2]
    
    # direct prediction
    DAR_nearest_neighbor <- R[,nni]*area[i]
    DAR_highest_rho <- R[,ri]*area[i]
    QPPQ_nearest_neighbor <- nn_quantile(U[,nni])*area[i]
    QPPQ_highest_rho <- nn_quantile(U[,ri])*area[i]
    
    # pseudo bivariate copula prediction
    corr_weighted_runoff_nn <- (sigma[i,nni]*R[,nni])*area[i]
    corr_weighted_runoff_rho <- (sigma[i,ri]*R[,ri])*area[i]
    
    # bivariate copula prediction
    bivar_norm_cop_nn <- nn_quantile(pnorm(sigma[i,nni]*Z[,nni]))*area[i]
    bivar_norm_cop_rho <- nn_quantile(pnorm(sigma[i,ri]*Z[,ri]))*area[i]
    # var_est <- 1-(sigma[i,ri]^2)
    # bivariate <- quantile(target,pnorm(mean_est))
    # bivar_low <- quantile(target, pnorm(mean_est - (2*sqrt(var_est))))
    # bivar_high <- quantile(target, pnorm(mean_est + (2*(sqrt(var_est)))))
    
    # inverse distance weighted
    w <- 1/(dmat[-i,i]^2)
    IDW_runoff <- ((R[,-i]%*%w)/sum(w))*area[i]
    IDW_log_runoff <- 10^((log10(R[,-i])%*%w)/sum(w))*area[i]
    IDW_nep <- nn_quantile((U[,-i]%*%w)/sum(w))*area[i]
    
    # pseudo multivariate copula prediction
    #sigma_R <- cov(R, method="spearman")
    sigma_R <- sigma
    rx <- R[,-i]
    mux <- apply(rx,2,mean)
    sigma_xx_R <- sigma_R[-i,-i]
    sigma_yx_R <- sigma_R[i,-i]
    inv_sigma_xx_R <- solve(sigma_xx_R)
    ry <- mean(R[,i]) + t((sigma_yx_R%*%inv_sigma_xx_R)%*%(t(rx[,])-mux))
    corr_weighted_runoff <- ry*area[i]
    
    # multivariate copula prediction 
    zx <- Z[,-i]
    sigma_xx <- sigma[-i,-i]
    sigma_yx <- sigma[i,-i]
    inv_sigma_xx <- solve(sigma_xx)
    zy <- t((sigma_yx%*%inv_sigma_xx)%*%t(zx))
    py <- pnorm(zy)
    multivar_norm_cop <- nn_quantile(pnorm(zy))*area[i]
    # var_y <- sigma[i,i]-(sigma_yx%*%inv_sigma_xx%*%t(sigma_yx)[,])
    # multivariate <- quantile(target,pnorm(zy))
    # multivar_low <- quantile(target, pnorm(zy - (2*sqrt(var_y)[[1]])))
    # multivar_high <- quantile(target, pnorm(zy + (2*sqrt(var_y)[[1]])))
    
    # mean of R, log(R), and U
    avg_runoff <- apply(R[,-i],1,mean)*area[i]
    avg_log_runoff <- 10^apply(log10(R[,-i]),1,mean)*area[i]
    avg_nep <- nn_quantile(apply(U[,-i],1,mean))*area[i]
    
    # combine estimates
    result <- data.frame(date = date,
                         target = sites[i],
                         donor_rho = sites[ri],
                         rho = sigma[i,ri],
                         observed = Q[,i],
                         DAR_nearest_neighbor,
                         DAR_highest_rho,
                         QPPQ_nearest_neighbor,
                         QPPQ_highest_rho,
                         corr_weighted_runoff_nn,
                         corr_weighted_runoff_rho,
                         bivar_norm_cop_nn,
                         bivar_norm_cop_rho,
                         IDW_runoff,
                         IDW_log_runoff,
                         IDW_nep,
                         corr_weighted_runoff,
                         multivar_norm_cop,
                         avg_runoff,
                         avg_log_runoff,
                         avg_nep,
                         Utarget = U[,i],
                         stringsAsFactors = F)
    
    result_all <- rbind(result_all,result)
    
  }
  
  rownames(result_all) <- NULL
  
  
  return(result_all)
}

# plot changes in NSE
plot_delta_nse <- function(estimates_best_e, estimates_best_w,estimates_fdc_e,estimates_fdc_w,estimates_cor_e,estimates_cor_w, estimates_fdc_cor_e,estimates_fdc_cor_w,models){
  
  estimates_best_e <- estimates_best_e %>%
    mutate(location="Mobile-Tombigbee",type="obs FDC + corr")
  
  estimates_best_w <- estimates_best_w %>%
    mutate(location="Galveston-Trinity",type="obs FDC + corr")
  
  all_est <- bind_rows(estimates_best_e, estimates_best_w,estimates_fdc_e,estimates_fdc_w,
                       estimates_cor_e,estimates_cor_w, estimates_fdc_cor_e,estimates_fdc_cor_w) %>%
    select(-c(date, donor_rho, rho, Utarget)) %>%
    gather(model,value,-observed,-target,-location,-type) %>%
    group_by(model,target,location,type) %>%
    summarize(NSE = nse(value,observed)) %>%
    ungroup() %>%
    mutate(model = forcats::fct_reorder(model,NSE),
           type = forcats::fct_relevel(type,c("est FDC + corr","est FDC","est corr","obs FDC + corr"))) %>%
    filter(model %in% models) 
  
  cols <- c("#c03728", "#919c4c","#828585", "#efe1c6","#f5c04a","dodgerblue4")
  
  ggplot(all_est) +
    geom_boxplot(aes(x=type,y=NSE,fill=model),outlier.size = 0.1) +
    facet_wrap(~location,ncol=1) +
    scale_fill_manual(values=cols,name="") +
    labs(x="") +
    theme_bw(base_size=14) +
    coord_cartesian(ylim=c(-1,1)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "top")
}

# plot time series for each method
plot_all_ts <- function(estimates_best_e, estimates_best_w,estimates_fdc_e,estimates_fdc_w,estimates_cor_e,estimates_cor_w,estimates_fdc_cor_e,estimates_fdc_cor_w,methods){
  
  estimates_best_e <- estimates_best_e %>%
    mutate(location="Mobile-Tombigbee",type="obs FDC + corr")
  
  estimates_best_w <- estimates_best_w %>%
    mutate(location="Galveston-Trinity",type="obs FDC + corr")
  
  all_ts <- bind_rows(estimates_best_e, estimates_best_w,estimates_fdc_e,estimates_fdc_w,
                       estimates_cor_e,estimates_cor_w, estimates_fdc_cor_e,estimates_fdc_cor_w) %>%
    filter(target %in% c("02464360","08066300") & date >= "2009-03-01" & date <="2009-09-01") %>%
    select(-c(donor_rho, rho, Utarget)) %>%
    gather(model,value,-target,-date,-observed,-location,-type) %>%
    mutate(type = forcats::fct_relevel(type,c("est FDC + corr","est FDC","est corr","obs FDC + corr"))) %>%
    filter(model %in% methods) 
  
  cols <- c("#c03728","#efe1c6","#919c4c","#f5c04a")
  
  ggplot(all_ts) +
    geom_line(aes(date,observed),color="grey30") +
    geom_line(aes(date,value,color=type)) +
    scale_color_manual(values=cols,name="") +
    labs(y=expression(m^3~s^{-1})) +
    facet_grid(model~location) +
    scale_y_log10() +
    theme_bw(base_size=14) +
    theme(legend.position = "top")
}
