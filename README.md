Data and Scripts for the paper "Empirical Bayesian interfrequency and interintensity correlations of horizontal Fourier amplitude spectra using NGA-West2 data" by Campbell,KW, Kuehn, N, and Bozorgnia, Y, accepted for publication in Earthquake Spectra.

data:
  - CB14_PSA_R0_R500_DATA.csv contains the data set used for development of the CB14 GMM
  - TotResidAllPer.csv contans total EAS residuals of the CB25 EAS GMM
  - TotResidAllPer_IM_ModES.csv contains total PSA and other IM residuals.
  - cor_ba19.csv contains correlation coefficients from BA19 of the target frequencies in this study with respect to 5Hz.

stan:
  - gmm_partition_corrre_cond_miss_tauM_phiM.stan contains stan code to calculate a bivariate mixed-effects model on total residuals with a uniform prior on the correlation coefficients.
  - gmm_partition_corrre_cond_miss_tauM_phiM_priorz.stan contains stan code to calculate a bivariate mixed-effects model on total residuals with a normal prior distrbution on the z-transformed correlation coefficients.

r:
  - EAS_EAS_Mdep_SDonly.R contains a script that calculates bivariate correlation coefficients for EAS total residuals.
  - stan_correlation_eas_priorsenstivity.R contains a script to perform prior sensitivity.
