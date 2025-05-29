################################################################################
# Script to per form prior sensitvity of correlation coefficients

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

inverse_ztransform <- function(rho) {
  return((exp(2 * rho) - 1) / (exp(2 * rho) + 1))
}
ztransform <- function(rho) {
  0.5 * log((1+rho) / (1-rho)) 
}

################################################################################

dir_base     <- "."
dir_data     <- file.path(dir_base, "../data/")
dir_stan     <- file.path(dir_base, "../stan/")
dir_ngaw2    <- file.path(dir_base, "../data/")


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
mod1 <- cmdstan_model(file.path(dir_stan, 'gmm_partition_corrre_cond_miss_tauM_phiM.stan'))
mod_z <- cmdstan_model(file.path(dir_stan, 'gmm_partition_corrre_cond_miss_tauM_phiM_priorz.stan'))

partition_data_rho_eas <- function(k1, k2,
                                   iter_warmup = 500, iter_sampling = 1000, chains = 4,
                                   model = 'base',
                                   pars_z_eq = c(0.55, 0.5),
                                   pars_z_stat = c(0.55, 0.5),
                                   pars_z_rec = c(0.55, 0.5),
                                   sd_z = 'given') {
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
  
  data_obs <- data_used[idx_obs_2,]
  if(sd_z != 'given') {
    pars_z_eq[2] <- 2/sqrt(length(unique(data_obs$EQID)) - 3)
    pars_z_stat[2] <- 2/sqrt(length(unique(data_obs$SSN)) - 3)
    pars_z_rec[2] <- 2/sqrt(nrow(data_obs) - 3)
  }
  
  mb <- c(4.5,5.5)
  m1_rec <- 1 * (data_used$M < mb[2]) - (data_used$M - mb[1]) / (mb[2] - mb[1]) * (data_used$M > mb[1] & data_used$M < mb[2])
  m2_rec <- 1 * (data_used$M >= mb[2]) + (data_used$M - mb[1]) / (mb[2] - mb[1]) * (data_used$M > mb[1] & data_used$M < mb[2])
  
  mageq <- unique(data_used[,c('M','EQID')])$M
  m1_eq <- 1 * (mageq < mb[2]) - (mageq - mb[1]) / (mb[2] - mb[1]) * (mageq > mb[1] & mageq < mb[2])
  m2_eq <- 1 * (mageq >= mb[2]) + (mageq - mb[1]) / (mb[2] - mb[1]) * (mageq > mb[1] & mageq < mb[2])
  
  data_list <- list(
    N = nrow(y_target),
    NEQ = max(eq),
    NSTAT = max(stat),
    N_obs_var2 = length(idx_obs_2),
    Y = y_target,
    eq = eq,
    stat = stat,
    idx_obs_var2 = idx_obs_2,
    M1_rec = m1_rec,
    M2_rec = m2_rec,
    M1_eq = m1_eq,
    M2_eq = m2_eq,
    pars_z_eq = pars_z_eq,
    pars_z_stat = pars_z_stat,
    pars_z_rec = pars_z_rec
  )
  
  if(model == 'base') {
    fit_partition <- mod1$sample(
      data = data_list,
      seed = 1701,
      parallel_chains = 3,
      chains = chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      refresh = 50,
      show_exceptions = FALSE
    )
  } else {
    fit_partition <- mod_z$sample(
      data = data_list,
      seed = 1701,
      parallel_chains = 3,
      chains = chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      refresh = 50,
      show_exceptions = FALSE
    )
  }
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
    subset(draws, variable = c('phi','tau','rho','z_'), regex = TRUE),
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
  
  return(list(rv = tmp, summary = tmp2, data = data_list))
}


# read correlations of BA19 model, all with respect to 5Hz
cor_ba19 <- read.csv(file.path(dir_data, 'cor_ba19.csv'))


results_cor_rv <- list()
results_cor_summary <- list()

results_cor_z_rv <- list()
results_cor_z_summary <- list()

results_cor_z2_rv <- list()
results_cor_z2_summary <- list()
results_cor_z2_data <- list()


n_chains <- 4
n_wp <- 500
n_samples <- 1000
n_post <- n_chains * n_wp

# not run
k1 = 14 # for 5Hz
k <- 1
for(k2 in 1:n_target) {
  if(k2 == k1) {next}
  
  print(c(names_target[k1], names_target[k2]))
  # non-informative prior for correlation coefficents
  tmp <- partition_data_rho_eas(k1, k2, iter_warmup = 500, iter_sampling = 1000, chains = 4, model = 'base')
  results_cor_rv[[k]] <- tmp$rv
  results_cor_summary[[k]] <- tmp$summary
  
  # weakly informative prior for correlation coefficients
  tmp <- partition_data_rho_eas(k1, k2, iter_warmup = 500, iter_sampling = 1000, chains = 4, model = 'z')
  results_cor_z_rv[[k]] <- tmp$rv
  results_cor_z_summary[[k]] <- tmp$summary
  
  # informative prior for correlation coefficients based on BA19
  z_eq <- 0.5 * log((1 + cor_ba19$eq[k2]) / (1 - cor_ba19$eq[k2]))
  z_stat <- 0.5 * log((1 + cor_ba19$stat[k2]) / (1 - cor_ba19$stat[k2]))
  z_rec <- 0.5 * log((1 + cor_ba19$rec[k2]) / (1 - cor_ba19$rec[k2]))
  
  tmp <- partition_data_rho_eas(k1, k2, iter_warmup = 500, iter_sampling = 1000, chains = 4, model = 'z',
                                pars_z_eq = c(z_eq,NA), pars_z_stat = c(z_stat,NA), pars_z_rec = c(z_rec,NA),
                                sd_z = 'calc')
  results_cor_z2_rv[[k]] <- tmp$rv
  results_cor_z2_summary[[k]] <- tmp$summary
  results_cor_z2_data[[k]] <- tmp$data
  
  k <- k + 1
}

# not run here
# save(results_cor_rv, results_cor_z_rv, results_cor_z2_rv,
#      results_cor_summary, results_cor_z_summary, results_cor_z2_summary,
#      results_cor_z2_data,
#      file = file.path(dir_result, "results_cor_priorsensitiviy.Rdata"))



length(results_cor_z2_data)
results_cor_z2_data[[idx1]]


make_plot_rho_post <- function(idx1) {
  
  p1 <- data.frame(unif_posterior = as_draws_matrix(results_cor_rv[[idx1]]$rho_eq),
                   prior1_posterior = as_draws_matrix(results_cor_z_rv[[idx1]]$rho_eq),
                   prior2_posterior = as_draws_matrix(results_cor_z2_rv[[idx1]]$rho_eq),
                   prior1_prior = tanh(rnorm(n_post, mean = 0.55, sd  = 0.5)),
                   prior2_prior = tanh(rnorm(n_post, mean = results_cor_z2_data[[idx1]]$pars_z_eq[1],
                                             sd  = results_cor_z2_data[[idx1]]$pars_z_eq[2]))
  ) %>% set_names(c('unif_posterior','prior1_posterior','prior2_posterior','prior1_prior','prior2_prior')) %>%
    pivot_longer(everything(), names_to = c('model','type'), names_sep = '_') %>%
    ggplot(aes(x = value, color = model, linetype = type)) +
    geom_density(linewidth = 1.5) +
    scale_color_manual(values = c('red','blue', 'black')) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', linewidth = 1.5) +
    labs(x = 'rho_dB')
  
  p2 <- data.frame(unif_posterior = as_draws_matrix(results_cor_rv[[idx1]]$rho_stat),
                   prior1_posterior = as_draws_matrix(results_cor_z_rv[[idx1]]$rho_stat),
                   prior2_posterior = as_draws_matrix(results_cor_z2_rv[[idx1]]$rho_stat),
                   prior1_prior = tanh(rnorm(n_post, mean = 0.55, sd  = 0.5)),
                   prior2_prior = tanh(rnorm(n_post, mean = results_cor_z2_data[[idx1]]$pars_z_stat[1],
                                             sd  = results_cor_z2_data[[idx1]]$pars_z_stat[2]))
  ) %>% set_names(c('unif_posterior','prior1_posterior','prior2_posterior','prior1_prior','prior2_prior')) %>%
    pivot_longer(everything(), names_to = c('model','type'), names_sep = '_') %>%
    ggplot(aes(x = value, color = model, linetype = type)) +
    geom_density(linewidth = 1.5) +
    scale_color_manual(values = c('red','blue', 'black')) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', linewidth = 1.5) +
    theme(legend.position = 'none') +
    labs(x = 'rho_dS')
  
  p3 <- data.frame(unif_posterior = as_draws_matrix(results_cor_rv[[idx1]]$rho_rec),
                   prior1_posterior = as_draws_matrix(results_cor_z_rv[[idx1]]$rho_rec),
                   prior2_posterior = as_draws_matrix(results_cor_z2_rv[[idx1]]$rho_rec),
                   prior1_prior = tanh(rnorm(n_post, mean = 0.55, sd  = 0.5)),
                   prior2_prior = tanh(rnorm(n_post, mean = results_cor_z2_data[[idx1]]$pars_z_rec[1],
                                             sd  = results_cor_z2_data[[idx1]]$pars_z_rec[2]))
  ) %>% set_names(c('unif_posterior','prior1_posterior','prior2_posterior','prior1_prior','prior2_prior')) %>%
    pivot_longer(everything(), names_to = c('model','type'), names_sep = '_') %>%
    ggplot(aes(x = value, color = model, linetype = type)) +
    geom_density(linewidth = 1.5) +
    scale_color_manual(values = c('red','blue', 'black')) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', linewidth = 1.5) +
    theme(legend.position = 'none') +
    labs(x = 'rho_dWS')
  
  p4 <- ggpubr::get_legend(p1)
  
  p1 <- p1 + theme(legend.position = 'none')
  
  pl <- p1 + patchwork::plot_spacer() + p2 + p3 + patchwork::plot_spacer() + p4 +
    patchwork::plot_layout(widths = c(1,0.01, 1)) +
    patchwork::plot_annotation(title = sprintf('Correlation between F1 = %.3fHz and F2 = %.3fHz',
                                               1/convert_to_numeric(results_cor_rv[[idx1]]$target1),
                                               1/convert_to_numeric(results_cor_rv[[idx1]]$target2)))
  
  return(pl)
}

idx1 <- 8 # 1 Hz
