library(purrr);library(dplyr)

list.files("sim_results/") %>% 
  map_df(~readRDS(paste0("sim_results/",.))) %>%
  group_by(s_true,weighting, scale, y_scale) %>% 
  summarise(bias = mean(m_01_bias),
            mu_v = mean(mu_v),
            mu_l = mean(mu_l),
            DoM = mean(DoM),
            ASMD = mean(ASMD)) %>% 
  saveRDS(paste0("combined_results/",format(Sys.Date(),"%m%d%Y"),"_calibration_sim_results.rds"))