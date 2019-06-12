### Setup
source("data_example/code/analysis_functions.R")

### Load data
data = readRDS("data_example/data/merged_data.rds")

### Alter the data so that BMI and race are not as different ??
# data$BMI[which(data$S == 0)] = data$BMI[which(data$S == 0)] + 2
# data$black[which(data$S == 0)] = rbinom(length(data$black[which(data$S == 0)]), 1, 0.20)

example_results = analysis(data,"6mo",trim_pctile = .9)$results %>% 
  left_join(analysis(data,"18mo",trim_pctile = .9)$results, by = "Method")
names(example_results)[2:7] = rep(c("$\\mu_0^{v}$","$\\mu_0^{rct}$","$\\text{bias}_{\\mu_0}$"),2)