
# load pacakges and functions
source("scripts/Worland/copulas/cop_utils.R")
pacman::p_load(tidyverse,stringr,copula,VineCopula,sf,plotly,geosphere,ggridges,mvnfast,ggsn,cowplot,forcats,Matrix, ggradar)

# load data
restore_hucs <- remake::fetch('restore_pp', remake_file='prepare_data.yml')
restore_hucs <- st_read("data/shapefiles/restore_hucs/restoreHUC12s.shp")
gage_all <- remake::fetch('gage_all', remake_file='prepare_data.yml') %>%
  rename(lat=dec_lat_va,lon=dec_long_va)
dv_list <- remake::fetch('dv_list', remake_file='predict_dv.yml')
est_fdcs <- remake::fetch('all_est_fdc', remake_file='predict_fdc.yml') 

# build huc4 basin dataset
east_huc4 <- east_data(gage_all)
west_huc4 <- west_data(gage_all)
east_dvs <- subset_dv(east_huc4, dv_list)
west_dvs <- subset_dv(west_huc4, dv_list)

# make map of study area
huc4_map(restore_hucs,east_huc4,west_huc4)

# plot time series of NEPs
plot_ts(east_dvs,west_dvs,gage_all,yr=2005)

# matrix of flow data
Q_e <- huc4_discharge(east_dvs)
Q_w <- huc4_discharge(west_dvs)

# vector of basin area
area_e <- east_huc4$basin_area
area_w <- west_huc4$basin_area

# matrix of unit discharge (runoff)
R_e <- t(t(Q_e)/area_e)
R_w <- t(t(Q_w)/area_w)

# covariance matrix of flow data
sigma_e <- cor(R_e, method="spearman")
sigma_w <- cor(R_w, method="spearman")

# estimated correlations
sigma_e_k <- read.table("data/farmer_coords/eastern_krige_rhos.csv", sep = ",") %>%
  set_names(paste0("0",parse_number(colnames(.)))) %>%
  as.matrix() %>%
  t()

diag(sigma_e_k) <- 1

sigma_w_k <- read.table("data/farmer_coords/western_krige_rhos.csv", sep = ",") %>%
  set_names(paste0("0",parse_number(colnames(.)))) %>%
  as.matrix() %>%
  t()

diag(sigma_w_k) <- 1

# matrix of NEPs for flow data
U_e <- apply(R_e,2,pobs)
U_w <- apply(R_w,2,pobs)

# matrix of z-scores for U
Z_e <- apply(U_e,2,qnorm)
Z_w <- apply(U_w,2,qnorm)

# list of coordinates for sites
coords_e <- st_coordinates(east_huc4)
coords_w <- st_coordinates(west_huc4)

# character vector of site numbers
sites_e <- unique(east_dvs$site_no)
sites_w <- unique(west_dvs$site_no)

# date vector
date <- unique(east_dvs$date)

# calculate best case scenario for all methods
estimates_best_e <- model_best(Q_e, R_e, U_e, Z_e, sigma_e, date, sites_e, coords_e, area_e)
estimates_best_w <- model_best(Q_w, R_w, U_w, Z_w, sigma_w, date, sites_w, coords_w, area_w)
plot_best_nse(estimates_best_e, estimates_best_w)
plot_ts_best_estimate(estimates_best_e,estimates_best_w)

# select methods for futher analysis
models <- c("multivar_norm_cop","bivar_norm_cop_rho","QPPQ_highest_rho",
             "corr_weighted_runoff","IDW_nep","IDW_log_runoff")

# compare log-likelihoods
log_lik_compare_e <- compare_log_lik(Z_e, sigma_e)
log_lik_compare_w <- compare_log_lik(Z_w, sigma_w)
log_lik_compare_e$df_plot
log_lik_compare_w$df_plot
df_e <- log_lik_compare_e$df
df_w <- log_lik_compare_w$df

# compare conditional variances
var_compare_e <- conditional_variance(Q_e, R_e, U_e, Z_e, sigma_e, date, sites_e, area_e, df=df_e)
var_compare_w <- conditional_variance(Q_w, R_w, U_w, Z_w, sigma_w, date, sites_w, area_w, df=df_w)
var_across_u(var_compare_w)
plot_ts_ci(var_compare_e)
check_percentiles(var_compare_e,var_compare_w)

# donor best ---> worst
rho_compare_e <- compare_across_rhos(Q_e, R_e, U_e, Z_e, sigma_e, date, sites_e, coords_e, area_e)
rho_compare_w <- compare_across_rhos(Q_w, R_w, U_w, Z_w, sigma_w, date, sites_w, coords_w, area_w)
plot_rho_compare(rho_compare_e, rho_compare_w, by=3)
plot_rho_compare2(rho_compare_e, rho_compare_w)

# fit copulas
cops_e <- fit_cops(U_e,sigma_e)
cops_w <- fit_cops(U_w,sigma_w)

# plot loglik
pll <- plot_logliks(cops_e,cops_w)

# check for tail dependence
tail_dep(U_e,U_w,p=0.01,pll)

# plot copula simulations
plot_sim_cops(N=2000,alpha=0.07)

# est fdcs and correlations
estimates_fdc_e <- model_fdc_est(Q_e, R_e, U_e, Z_e, sigma_e, date, sites_e, 
                                 coords_e, area_e, est_fdcs) %>%
  mutate(location="Mobile-Tombigbee",type="est FDC")

estimates_fdc_w <- model_fdc_est(Q_w, R_w, U_w, Z_w, sigma_w, date, sites_w, 
                                 coords_w, area_w, est_fdcs) %>%
  mutate(location="Galveston-Trinity",type="est FDC")

estimates_cor_e <- model_fdc_est(Q_e, R_e, U_e, Z_e, sigma_e_k, date, sites_e, 
                              coords_e, area_e,est_fdcs) %>%
  mutate(location="Mobile-Tombigbee",type="est corr")

estimates_cor_w <- model_fdc_est(Q_w, R_w, U_w, Z_w, sigma_w_k, date, sites_w, 
                              coords_w, area_w,est_fdcs) %>%
  mutate(location="Galveston-Trinity",type="est corr")

estimates_fdc_cor_e <- model_fdc_est(Q_e, R_e, U_e, Z_e, sigma_e_k, date, sites_e, 
                                     coords_e, area_e, est_fdcs) %>%
  mutate(location="Mobile-Tombigbee",type="est FDC + corr")

estimates_fdc_cor_w <- model_fdc_est(Q_w, R_w, U_w, Z_w, sigma_w_k, date, sites_w, 
                                     coords_w, area_w, est_fdcs) %>%
  mutate(location="Galveston-Trinity",type="est FDC + corr")

# boxlplot of all methods, location, and type
plot_delta_nse(estimates_best_e, estimates_best_w,estimates_fdc_e,estimates_fdc_w,estimates_cor_e,
               estimates_cor_w,estimates_fdc_cor_e,estimates_fdc_cor_w,models)

plot_all_ts(estimates_best_e, estimates_best_w,estimates_fdc_e,estimates_fdc_w,estimates_cor_e,
            estimates_cor_w,estimates_fdc_cor_e,estimates_fdc_cor_w,models)

  
# t(gage_all[which(gage_all$site_no=="02430880" & gage_all$decade==2000),11:37])*35.3147
# t(gage_all[which(gage_all$site_no=="02430880" & gage_all$decade==2000),11:37])