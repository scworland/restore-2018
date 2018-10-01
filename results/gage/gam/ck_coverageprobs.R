library(feather)


L1df <- read_feather("all_gage_looest_L1.feather")
C <- L1df[L1df$in_model_L1 == 1,]
n <- length(C$site_no)
m <-    sum(C$L1 < C$est_upr_L1     & C$L1 > C$est_lwr_L1)
mloo <- sum(C$L1 < C$loo_est_upr_L1*bc & C$L1 > C$loo_est_lwr_L1*(1/bc))
message("             L1: coverage probability=", sprintf("%.3f", round(m/n,    digits=6)), " and ",
                     "LOO coverage probability=", sprintf("%.3f", round(mloo/n, digits=6)))

T2df <- read_feather("all_gage_looest_T2.feather")
C <- T2df[T2df$in_model_T2 == 1,]
n <- length(C$site_no)
m <-    sum(C$T2 < C$est_upr_T2     & C$T2 > C$est_lwr_T2)
mloo <- sum(C$T2 < C$loo_est_upr_T2*bc & C$T2 > C$loo_est_lwr_T2*(1/bc))
message("             T2: coverage probability=", sprintf("%.3f", round(m/n,    digits=6)), " and ",
                     "LOO coverage probability=", sprintf("%.3f", round(mloo/n, digits=6)))

T3df <- read_feather("all_gage_looest_T3.feather")
C <- T3df[T3df$in_model_T3 == 1,]
n <- length(C$site_no)
m <-    sum(C$T3 < C$est_upr_T3     & C$T3 > C$est_lwr_T3)
mloo <- sum(C$T3 < C$loo_est_upr_T3*bc & C$T3 > C$loo_est_lwr_T3*(1/bc))
message("             T3: coverage probability=", sprintf("%.3f", round(m/n,    digits=6)), " and ",
                     "LOO coverage probability=", sprintf("%.3f", round(mloo/n, digits=6)))

T4df <- read_feather("all_gage_looest_T4.feather")
C <- T4df[T4df$in_model_T4 == 1,]
n <- length(C$site_no)
m <-    sum(C$T4 < C$est_upr_T4     & C$T4 > C$est_lwr_T4)
mloo <- sum(C$T4 < C$loo_est_upr_T4*bc & C$T4 > C$loo_est_lwr_T4*(1/bc))
message("             T4: coverage probability=", sprintf("%.3f", round(m/n,    digits=6)), " and ",
                     "LOO coverage probability=", sprintf("%.3f", round(mloo/n, digits=6)))

T5df <- read_feather("all_gage_looest_T5.feather")
C <- T5df[T5df$in_model_T5 == 1,]
n <- length(C$site_no)
m <-    sum(C$T5 < C$est_upr_T5     & C$T5 > C$est_lwr_T5)
mloo <- sum(C$T5 < C$loo_est_upr_T5*bc & C$T5 > C$loo_est_lwr_T5*(1/bc))
message("             T5: coverage probability=", sprintf("%.3f", round(m/n,    digits=6)), " and ",
                     "LOO coverage probability=", sprintf("%.3f", round(mloo/n, digits=6)))

bc <- 1
#bc <- 1.01 # a preliminary screening shows bulk off of the LOO by about 1%.
PPLOdf <- read_feather("all_gage_looest_pplo.feather")
C <- PPLOdf[PPLOdf$in_model_pplo == 1,]
n <- length(C$site_no)
m <-    sum(C$flowtime < C$est_upr_flowtime     & C$flowtime > C$est_lwr_flowtime)
mloo <- sum(C$flowtime < C$loo_est_upr_flowtime*bc & C$flowtime > C$loo_est_lwr_flowtime*(1/bc))
message("PPLO (flowtime): coverage percentage=", sprintf("%.3f", round(m/n,    digits=4)), " and ",
                     "LOO coverage percentage=", sprintf("%.3f", round(mloo/n, digits=4)))

