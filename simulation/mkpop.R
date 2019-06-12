library(dplyr);library(purrr);library(MASS);library(stringr)
source("sim_functions.R")
source("utils.R")
source("config.R")

start_time = Sys.time()

### Set the seed
set.seed(13783)

### Make the population data files
if((!"pop_data_params.rds" %in% list.files("pop_data/")) |
   (!identical(pop_data_params, readRDS("pop_data/pop_data_params.rds")))){
  
  ### Make "x" matrix:
  x = mvrnorm(1000000, c(0,0,0,0), Sigma = matrix(c(1,0,0,0,
                                                    0,1,0,0,
                                                    0,0,1,0,
                                                    0,0,0,1),nrow=4))
  
  seq_along(pop_data_params) %>% 
    map(function(i){
      pop = data_pop(x, pop_data_params[[i]]$scale, 1, 0, 0.5, 2, pop_data_params[[i]]$s_true)
      saveRDS(pop, paste0("pop_data/pop_",i,".rds"))
    })
  
  saveRDS(pop_data_params,"pop_data/pop_data_params.rds")
} else{
  cat("Population data in config are already made!")
}

end_time = Sys.time()
print(difftime(end_time,start_time,units="mins"))