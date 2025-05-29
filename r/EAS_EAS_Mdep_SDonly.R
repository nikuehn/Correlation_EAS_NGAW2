##*************************************************************************##
##***** BIVARIATE BAYESIAN INTERFREQUENCY CORRELATION ANALYSIS OF EAS *****##
##*************************************************************************##
##*****      MAGNITUDE-INDEPENDENT COMPONENT STANDARD DEVIATIONS      *****##
##*************************************************************************##
##*****         MAGNITUDE-DEPENDENT TOTAL STANDARD DEVIATIONS         *****##
##*************************************************************************##
##*****                MODEL ES - EVENT AND SITE TERMS                *****##
##*************************************************************************##

IMtype1 = "EAS"
IMtype2 = "EAS"
Mod     = "ES"

library(cmdstanr)
library(posterior)
library(bayesplot)
library(matrixStats)
library(tidyverse)


################################################################################
`%notin%` <- Negate(`%in%`)

convert_to_numeric <- function(str) {
  str <- sub("^T", "", str)
  str <- sub("p", ".", str)
  as.numeric(str)
}

################################################################################

dir_base     <- "."
dir_data     <- file.path(dir_base, "../data/")
dir_stan     <- file.path(dir_base, "../stan/")
dir_ngaw2    <- file.path(dir_base, "../data/")
dir_res      <- file.path(dir_base, "results", paste0("Mod",Mod),IMtype1)

################################################################################
# read NGA W2 flatfile and total FAS residuals
# total FAS residual file does not include SSN, no need to gt from NGA W2
ff_ngaw2 <- read.csv(paste0(dir_ngaw2, 'CB14_PSA_R0_R500_DATA.csv'))
totres_fas <- read.csv(paste0(dir_data, 'TotResidAllPer.csv'))
totres_fas_combined <- left_join(totres_fas,
                                 ff_ngaw2 %>%
                                   select(RSN, EQID, SSN) %>%
                                   set_names(c('RSN','EQID','SSN')),
                                 by = 'RSN')
dim(totres_fas_combined)

names_target <- names(totres_fas)[!is.na(str_extract(names(totres_fas),pattern = "T[0-9]"))]
per_target <- Round(convert_to_numeric(names_target), 3)
freqs_target <- Round(1/convert_to_numeric(names_target),3)
n_target <- length(names_target)

################################################################################
# define Stan model
mod <- cmdstan_model(file.path(dir_stan, 'gmm_partition_corrre_cond_miss_tauM_phiM.stan'))

partition_data_rho_eas <- function(k1, k2,
                                   iter_warmup = 500, iter_sampling = 1000, chains = 4) {
  data_sel1 <- totres_fas_combined %>%
    select(RSN, EQID, SSN, M, names_target[k1], names_target[k2])
  
  # find target variable which has larger number of data
  n1 <- sum(!is.na(data_sel1[,names_target[k1]]))
  n2 <- sum(!is.na(data_sel1[,names_target[k2]]))
  
  if(n1 > n2) {
    i1 <- 1
    i2 <- 2
    name_1 <- names_target[k1]
    name_2 <- names_target[k2]
  } else {
    i1 <- 2
    i2 <- 1
    name_1 <- names_target[k2]
    name_2 <- names_target[k1]
  }
  data_used <- data_sel1[!is.na(data_sel1[,name_1]),]
  print(dim(data_used))
  
  y_target <- data_used[,c(name_1, name_2)]
  idx_miss_2 <- which(is.na(y_target[,2]))
  idx_obs_2 <- which(!is.na(y_target[,2]))
  
  eq <- as.numeric(factor(data_used$EQID, levels = unique(data_used$EQID)))
  stat <- as.numeric(factor(data_used$SSN, levels = unique(data_used$SSN)))
  
  y_target[idx_miss_2,2] <- -999
  
  mb <- c(4.5,5.5)
  m1_rec <- 1 * (data_used$M < mb[2]) - (data_used$M - mb[1]) / (mb[2] - mb[1]) * (data_used$M > mb[1] & data_used$M < mb[2])
  m2_rec <- 1 * (data_used$M >= mb[2]) + (data_used$M - mb[1]) / (mb[2] - mb[1]) * (data_used$M > mb[1] & data_used$M < mb[2])
  
  mageq <- unique(data_used[,c('M','EQID')])$M
  m1_eq <- 1 * (mageq < mb[2]) - (mageq - mb[1]) / (mb[2] - mb[1]) * (mageq > mb[1] & mageq < mb[2])
  m2_eq <- 1 * (mageq >= mb[2]) + (mageq - mb[1]) / (mb[2] - mb[1]) * (mageq > mb[1] & mageq < mb[2])
  
  bin_eq <- cut(mageq, breaks = c(3,4.5,5.5,9),labels = FALSE)
  bin_rec <- cut(data_used$M, breaks = c(3,4.5,5.5,9),labels = FALSE)
  bin_labels <- levels(cut(c(3.5,5,7), breaks = c(3,4.5,5.5,9)))
  
  data_list <- list(
    N = nrow(y_target),
    NEQ = max(eq),
    NSTAT = max(stat),
    N_obs_var2 = length(idx_obs_2),
    N_bin = length(bin_labels),
    Y = y_target,
    eq = eq,
    stat = stat,
    idx_obs_var2 = idx_obs_2,
    M1_rec = m1_rec,
    M2_rec = m2_rec,
    M1_eq = m1_eq,
    M2_eq = m2_eq,
    bin_eq = bin_eq,
    bin_rec = bin_rec
  )
  
  fit_partition <- mod$sample(
    data = data_list,
    seed = 1701,
    parallel_chains = 4,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = 50,
    show_exceptions = FALSE
  )
  print(fit_partition$cmdstan_diagnose())
  print(fit_partition$diagnostic_summary())
  draws <- fit_partition$draws()
  rv <- as_draws_rvars(draws)
  
  rv$rho_total_sm <- (rv$phi_ss_sm[1] * rv$phi_ss_sm[2] * rv$rho_rec +
                        rv$phi_s2s[1] * rv$phi_s2s[2] * rv$rho_stat +
                        rv$tau_sm[1] * rv$tau_sm[2] * rv$rho_eq) /
    (sqrt(rv$phi_ss_sm[1]^2 + rv$phi_s2s[1]^2 + rv$tau_sm[1]^2) *
       sqrt(rv$phi_ss_sm[2]^2 + rv$phi_s2s[2]^2 + rv$tau_sm[2]^2))
  
  rv$rho_total_lm <- (rv$phi_ss_lm[1] * rv$phi_ss_lm[2] * rv$rho_rec +
                        rv$phi_s2s[1] * rv$phi_s2s[2] * rv$rho_stat +
                        rv$tau_lm[1] * rv$tau_lm[2] * rv$rho_eq) /
    (sqrt(rv$phi_ss_lm[1]^2 + rv$phi_s2s[1]^2 + rv$tau_lm[1]^2) *
       sqrt(rv$phi_ss_lm[2]^2 + rv$phi_s2s[2]^2 + rv$tau_lm[2]^2))
  
  tmp <- data.frame(target1 = name_1, 
                    target2 = name_2, 
                    rho_eq = rv$rho_eq,
                    rho_stat = rv$rho_stat, 
                    rho_rec = rv$rho_rec,
                    rho_total_sm = rv$rho_total_sm, rho_total_lm = rv$rho_total_lm,
                    tau_target1_sm = rv$tau_sm[1], tau_target2_sm = rv$tau_sm[2],
                    tau_target1_lm = rv$tau_lm[1], tau_target2_lm = rv$tau_lm[2],
                    phi_s2s_target1 = rv$phi_s2s[1], phi_s2s_target2 = rv$phi_s2s[2],
                    phi_ss_target1_sm = rv$phi_ss_sm[1], phi_ss_target2_sm = rv$phi_ss_sm[2],
                    phi_ss_target1_lm = rv$phi_ss_lm[1], phi_ss_target2_lm = rv$phi_ss_lm[2])
  
  tmp2 <- summarise_draws(
    subset(draws, variable = c('phi','tau','rho'), regex = TRUE),
    "mean",
    "median",
    "sd",
    "mad",
    "rhat",
    "ess_bulk",
    "ess_tail",
    ~quantile(.x, probs = c(0.025, 0.975))
  ) %>%
    set_names("variable","mean","median","sd","mad","rhat","ess_bulk","ess_tail","q2.5","q97.5")
  
  return(list(rv = tmp, summary = tmp2))
}


results_cor_rv <- list()
results_cor_summary <- list()
k <- 1
for(k1 in 1:n_target) {
  for(k2 in (k1+1):n_target) {
    if(k2 > length(names_target)) {break}
    
    print(c(names_target[k1], names_target[k2]))
    tmp <- partition_data_rho_eas(k1, k2, iter_warmup = 500, iter_sampling = 1000, chains = 4)
    results_cor_rv[[k]] <- tmp$rv
    results_cor_summary[[k]] <- tmp$summary
    k <- k + 1
  }
}
save(results_cor_rv, results_cor_summary, file=file.path(dir_res, 'results_',IMtype1,'_',IMtype2,'_Mdep_SDonly.Rdata'))

load(file=paste0(dir_res, 'results_',IMtype1,'_',IMtype2,'_Mdep_SDonly.Rdata'))
for(k in 1:length(results_cor_summary)) {
  results_cor_summary[[k]]$target1 <- results_cor_rv[[k]]$target1
  results_cor_summary[[k]]$target2 <- results_cor_rv[[k]]$target2
  results_cor_summary[[k]]$freq1 <- 1/convert_to_numeric(results_cor_rv[[k]]$target1)
  results_cor_summary[[k]]$freq2 <- 1/convert_to_numeric(results_cor_rv[[k]]$target2)
}

res_summary <- results_cor_summary[[1]]
for(k in 2:length(results_cor_summary)) {
  res_summary <- rbind(res_summary, results_cor_summary[[k]])
}
write.csv(res_summary, file=file.path(dir_res, "summary_",IMtype1,"_",IMtype2,"_Mdep_SDonly.csv"),
          row.names = FALSE)



func <- function(x) {
  c(mean(x), quantile(x, probs = c(0.025, 0.975)))
}
res_summary2 <- matrix(nrow = length(results_cor_rv), ncol = 3 * length(results_cor_rv[[1]]) - 4)
for(k in 1:length(results_cor_rv)) {
  res_summary2[k,] <- c(convert_to_numeric(results_cor_rv[[k]]$target1),
                        convert_to_numeric(results_cor_rv[[k]]$target2),
                        unlist(lapply(results_cor_rv[[k]][3:length(results_cor_rv[[k]])], func)))
}
names_summary2 <- unlist(lapply(names(results_cor_rv[[k]][3:length(results_cor_rv[[k]])]), function(x) {
  c(paste(x, 'mean', sep = '_'),
    paste(x, 'q025', sep = '_'),
    paste(x, 'q975', sep = '_'))
}))
df_summary2 <- data.frame(res_summary2) %>% set_names(c('T1','T2',names_summary2))
write.csv(df_summary2, file=file.path(dir_res,"results_",IMtype1,"_",IMtype2,"_Mdep_SDonly.csv"),
          row.names = names_target)

