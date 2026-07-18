# packages
library(dplyr)
library(jomo)
library(mitools)
library(mice)
library(foreach)
library(iterators)
library(parallel)
library(LaplacesDemon)

# get array ID as a string
index <- Sys.getenv("SLURM_ARRAY_TASK_ID")

# set random seed
seed <- 777 + as.numeric(index)
set.seed(seed = seed)

# get the command line to pass into script
args = commandArgs(trailingOnly=TRUE)

pr_cmiss <- as.numeric(args[1])
pr_e <- as.numeric(args[2])
pr_c <- as.numeric(args[3])
pr_y <- as.numeric(args[4])

rr_e_giv_c <- as.numeric(args[5])
rr_e <- as.numeric(args[6])
rr_c <- as.numeric(args[7])
pr_e_giv_c <- as.numeric(args[8])
pry_base <- as.numeric(args[9])

nsims = 2000
n_size = 2000

ncores <- as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))

ret1 <- parallel::mclapply(1:nsims, function(sim){
  
  print(sim)
  # data generation
  M <- rbinom(n = n_size, size = 1, prob = pr_cmiss)
  C <- rbinom(n = n_size, size = 1, prob = pr_c)
  temp_df <- data.frame(M = M, C = C)
  for (i in 1:nrow(temp_df)) {
    temp_df$E[i] <- rbinom(n = 1, size = 1, prob = rr_e_giv_c^temp_df$C[i]*pr_e_giv_c)
    temp_df$Y[i] <- rbinom(n = 1, size = 1, prob = rr_e^temp_df$E[i]*rr_c^temp_df$C[i]*pry_base)
  }
  
  # no missing method
  fit_nm <- tryCatch(glm(Y ~ E + C, data = temp_df, 
                         family = binomial(link = "log")), 
                     error = function(e) NA)
  
  if (sum(is.na(fit_nm)) == 1) {
    # break
    skip_setting_nm = TRUE

      beta_nm = NA
      beta_nm_se = NA
      lower_ci_nm = NA
      upper_ci_nm = NA
      capture_nm = NA
  }else{
  skip_setting_nm = FALSE
  beta_nm <- unname(fit_nm$coefficients[2])
  beta_nm_se <- coef(summary(fit_nm))[2, 2]
  lower_ci_nm <- coef(summary(fit_nm))[2, 1] - 1.96 * coef(summary(fit_nm))[2, 2]
  upper_ci_nm <- coef(summary(fit_nm))[2, 1] + 1.96 * coef(summary(fit_nm))[2, 2]
  capture_nm <- rr_e >= exp(lower_ci_nm) & rr_e <= exp(upper_ci_nm)
  }
  
  # mcim method
  temp_df$C_miss <- (1 - M)*C
  fit_mcim <- tryCatch(glm(Y ~ E + C_miss + M, data = temp_df, 
                           family = binomial(link = "log")), 
                       error = function(e) NA)
  
  # if fail to converge, skip this scenario
  if (sum(is.na(fit_mcim)) == 1) {
    # break
    skip_setting_mcim = TRUE
      beta_mcim = NA
      beta_mcim_se = NA
      lower_ci_mcim = NA
      upper_ci_mcim = NA
      capture_mcim = NA 
  }else{
  skip_setting_mcim = FALSE
  beta_mcim <- unname(fit_mcim$coefficients[2])
  beta_mcim_se <- coef(summary(fit_mcim))[2, 2]
  lower_ci_mcim <- coef(summary(fit_mcim))[2, 1] - 1.96 * coef(summary(fit_mcim))[2, 2]
  upper_ci_mcim <- coef(summary(fit_mcim))[2, 1] + 1.96 * coef(summary(fit_mcim))[2, 2]
  capture_mcim <- rr_e >= exp(lower_ci_mcim) & rr_e <= exp(upper_ci_mcim)
  }
  
  
  # multiple imputation method
  temp_df$C_mimp <- ifelse(temp_df$M == 1, 
                           yes = NA, no = C)
  temp_df_mimp <- temp_df %>% 
    dplyr::select(Y, E, C_mimp, M)
  temp_df_mimp <- data.frame(temp_df_mimp)
  
  out_temp_df <- jomo1(temp_df_mimp[,-c(4)], nimp = 10)
  out_temp_df <- subset(out_temp_df, Imputation>0)
  mi_list <- imputationList(split(out_temp_df, out_temp_df$Imputation))
  fit_mimp <- with(mi_list, tryCatch(glm(Y ~ E + C_mimp , 
                                         family = binomial(link = "log")), 
                                     error = function(e){
                                       NA
                                     }))
  
  
  # if fail to converge, skip this scenario
  if (sum(is.na(fit_mimp)) > 5) {
    # break
    skip_setting_mimp = TRUE
      beta_mimp = NA
      beta_mimp_se = NA
      lower_ci_mimp = NA
      upper_ci_mimp = NA
      capture_mimp = NA
  }else{
    if(length(which(is.na(fit_mimp)))>0){
      fit_mimp=fit_mimp[-which(is.na(fit_mimp))]
    }
    skip_setting_mimp = FALSE
    beta_mimp <- summary(pool(as.mira(fit_mimp)))[2,2]
    beta_mimp_se <- summary(pool(as.mira(fit_mimp)))[2,3]
    lower_ci_mimp <- summary(pool(as.mira(fit_mimp)))[2,2] - 1.96 * summary(pool(as.mira(fit_mimp)))[2,3]
    upper_ci_mimp <- summary(pool(as.mira(fit_mimp)))[2,2] + 1.96 * summary(pool(as.mira(fit_mimp)))[2,3]
    capture_mimp <- rr_e >= exp(lower_ci_mimp) & rr_e <= exp(upper_ci_mimp)
  }
  

  
  # complete case method
  temp_df2 <- temp_df[temp_df$M == 0,]
  
  # in extreme cases, skip this scenario if all Y's are 0 or 1
  if ((sum(temp_df2$Y) == nrow(temp_df2)) | (sum(temp_df2$E) == nrow(temp_df2)) | 
      (sum(temp_df2$C) == nrow(temp_df2))) {
    # break
    skip_setting_cc = TRUE
      beta_cc = NA
      beta_cc_se = NA
      lower_ci_cc = NA
      upper_ci_cc = NA
      capture_cc = NA
   
  }else{
  
    fit_cc <- tryCatch(glm(Y ~ E + C, data = temp_df2,
                         family = binomial(link = "log")),
                     error = function(f) NA)
  
  # if fail to converge, skip this scenario
    if (sum(is.na(fit_cc)) == 1) {
    # break
      skip_setting_cc = TRUE
       beta_cc = NA
      beta_cc_se = NA
      lower_ci_cc = NA
      upper_ci_cc = NA
      capture_cc = NA
  } else{
    skip_setting_cc = FALSE
    beta_cc <- unname(fit_cc$coefficients[2])
    beta_cc_se <- coef(summary(fit_cc))[2, 2]
    lower_ci_cc <- coef(summary(fit_cc))[2, 1] - 1.96 * coef(summary(fit_cc))[2, 2]
    upper_ci_cc <- coef(summary(fit_cc))[2, 1] + 1.96 * coef(summary(fit_cc))[2, 2]
    capture_cc <-  rr_e >= exp(lower_ci_cc) & rr_e <= exp(upper_ci_cc)
  }
 }

  
  return(ret1 <- list(
    pr_cmiss = pr_cmiss,
    pr_e = pr_e,
    pr_c = pr_c,
    pr_y = pr_y,
    rr_e_giv_c = rr_e_giv_c,
    rr_e = rr_e,
    rr_c = rr_c,
    skip_setting_nm = skip_setting_nm,
    skip_setting_mcim = skip_setting_mcim,
    skip_setting_mimp = skip_setting_mimp,
    skip_setting_cc = skip_setting_cc,
    beta_nm = beta_nm,
    beta_nm_se = beta_nm_se,
    lower_ci_nm = lower_ci_nm,
    upper_ci_nm = upper_ci_nm,
    capture_nm = capture_nm,
    beta_mcim = beta_mcim,
    beta_mcim_se = beta_mcim_se,
    lower_ci_mcim = lower_ci_mcim,
    upper_ci_mcim = upper_ci_mcim,
    capture_mcim = capture_mcim, 
    beta_mimp = beta_mimp,
    beta_mimp_se = beta_mimp_se,
    lower_ci_mimp = lower_ci_mimp,
    upper_ci_mimp = upper_ci_mimp,
    capture_mimp = capture_mimp,
    beta_cc = beta_cc,
    beta_cc_se = beta_cc_se,
    lower_ci_cc = lower_ci_cc,
    upper_ci_cc = upper_ci_cc,
    capture_cc = capture_cc
  ))
}, 
mc.cores = ncores, mc.silent = FALSE)

sim_result <- data.frame(do.call(rbind, ret1))
sim_result <- apply(sim_result, MARGIN = 2, function(x){
  unlist(x)
})
sim_result <- as.data.frame(sim_result)

# # what we want for each scenario
# sim_result2 <- data.frame(
#   pr_cmiss, pr_e, pr_c, pr_y, rr_e_giv_c, rr_e, rr_c,
#   beta_mean_nm = mean(sim_result$beta_nm, na.rm = T),
#   beta_med_nm = median(sim_result$beta_nm, na.rm = T),
#   beta_iqr_nm = IQR(sim_result$beta_nm, na.rm = T),
#   se_mean_nm = mean(sim_result$beta_nm_se, na.rm = T),
#   se_med_nm = median(sim_result$beta_nm_se, na.rm = T),
#   se_iqr_nm = IQR(sim_result$beta_nm_se, na.rm = T),
#   captured_in_ci_nm = mean(sim_result$capture_nm, na.rm = T),
#   mean_lower_ci_nm = mean(sim_result$lower_ci_nm, na.rm = T),
#   mean_upper_ci_nm = mean(sim_result$upper_ci_nm, na.rm = T),
#   beta_mean_mcim = mean(sim_result$beta_mcim, na.rm = T),
#   beta_med_mcim = median(sim_result$beta_mcim, na.rm = T),
#   beta_iqr_mcim = IQR(sim_result$beta_mcim, na.rm = T),
#   se_mean_mcim = mean(sim_result$beta_mcim_se, na.rm = T),
#   se_med_mcim = median(sim_result$beta_mcim_se, na.rm = T),
#   se_iqr_mcim = IQR(sim_result$beta_mcim_se, na.rm = T),
#   captured_in_ci_mcim = mean(sim_result$capture_mcim, na.rm = T),
#   mean_lower_ci_mcim = mean(sim_result$lower_ci_mcim, na.rm = T),
#   mean_upper_ci_mcim = mean(sim_result$upper_ci_mcim, na.rm = T),
#   beta_mean_mimp = mean(sim_result$beta_mimp, na.rm = T),
#   beta_med_mimp = median(sim_result$beta_mimp, na.rm = T),
#   beta_iqr_mimp = IQR(sim_result$beta_mimp, na.rm = T),
#   se_mean_mimp = mean(sim_result$beta_mimp_se, na.rm = T),
#   se_med_mimp = median(sim_result$beta_mimp_se, na.rm = T),
#   se_iqr_mimp = IQR(sim_result$beta_mimp_se, na.rm = T),
#   captured_in_ci_mimp = mean(sim_result$capture_mimp, na.rm = T),
#   mean_lower_ci_mimp = mean(sim_result$lower_ci_mimp, na.rm = T),
#   mean_upper_ci_mimp = mean(sim_result$upper_ci_mimp, na.rm = T),
#   beta_mean_cc = mean(sim_result$beta_cc, na.rm = T),
#   beta_med_cc = median(sim_result$beta_cc, na.rm = T),
#   beta_iqr_cc = IQR(sim_result$beta_cc, na.rm = T),
#   se_mean_cc = mean(sim_result$beta_cc_se, na.rm = T),
#   se_med_cc = median(sim_result$beta_cc_se, na.rm = T),
#   se_iqr_cc = IQR(sim_result$beta_cc_se, na.rm = T),
#   captured_in_ci_cc = mean(sim_result$capture_cc, na.rm = T),
#   mean_lower_ci_cc = mean(sim_result$lower_ci_cc, na.rm = T),
#   mean_upper_ci_cc = mean(sim_result$upper_ci_cc, na.rm = T)
# )


# output files
write.csv(sim_result, file = paste0("simulation_for_",index,".csv"))



