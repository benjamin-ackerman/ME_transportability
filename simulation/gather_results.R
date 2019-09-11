library(purrr);library(dplyr)
source("utils.R")

list.files("sim_results/") %>% 
  map_df(~readRDS(paste0("sim_results/",.))) %T>%
  saveRDS(., paste0("combined_results/",format(Sys.Date(),"%m%d%Y"),"_calibration_sim_results_full.rds")) %>% 
  group_by(s_true,weighting, scale, y_scale) %>% 
  summarise(bias = mean(m_01_bias),
            coverage = mean(coverage(mu_v, se_mu_v, mean(mu_l))),
            mu_v = mean(mu_v),
            mu_l = mean(mu_l),
            DoM = mean(DoM),
            ASMD = mean(ASMD),
            mean_pct_wt_outliers = mean(pct_wt_outliers)) %>% 
  saveRDS(paste0("combined_results/",format(Sys.Date(),"%m%d%Y"),"_calibration_sim_results.rds"))