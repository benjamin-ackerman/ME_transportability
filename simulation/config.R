library(dplyr);library(purrr)

################################
###### Parameters to Vary ######
################################

# Number of samples to draw from each population
n_samples = 1

# Parameters for population data: scale and true S models
pop_data_params = list(scale = seq(0, 1, by = .2), # VARY THIS
                       b1 = 1, b2 = 0, b3 = 0.5, b4 = 2, # relative S model coefficients
                       s_true = list(list(), # VARY THIS
                                     list("x3^2"),
                                     list("x4^2"),
                                     list("x3*x4"),
                                     list("x3^2","x4^2","x3*x4"),
                                     list("x1*x4"),
                                     list("x1*x3"))
                       ) %>% cross()

# Parameters for sample data: y_scale
sample_data_params = list(y_scale = seq(0, 1, by = 0.2), # VARY THIS
                          n = 1000, # Sample size of each validation data/trial
                          g0 = 0, g1 = 2, g2 = 0, sigma2_z = 1, # Z model
                          a1 = 1, a2 = 0, a3 = 0, a4 = 1, a5 = 2, a6 = 0.5 #Y model coefficients
                          ) %>% cross()
