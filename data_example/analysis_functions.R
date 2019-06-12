### Setup
library(dplyr);library(ggplot2);library(WeightIt);library(cobalt)
#source("simulation/code/utils.R")
source("utils.R")

fit_selection_model = function(data, method){
  df = data %>% 
    tidyr::drop_na(ID:education) 
  
  if(method == "gbm"){
    set.seed(13783)
    gbm.ATT <- weightit(S ~ BMI + as.factor(age_cat) + as.factor(sex) + as.factor(black) + as.factor(education),
                        data=df,
                        method="twang",
                        stop.method="ks.mean",
                        estimand="ATT")
    
    df$weight = gbm.ATT$weights
    df$weight[which(df$S == 1)] = 0
    df$ps = gbm.ATT$ps
  }
  
  if(method == "lr"){
    model = glm(S ~ BMI + as.factor(sex) + as.factor(black) + as.factor(age_cat) + as.factor(education), data = df, family = 'binomial')
    
    df$ps = predict(model, type = 'response')
    df$weight = ifelse(df$S == 0, df$ps/(1-df$ps), 0)
    
    gbm.ATT = NULL
  }
  
  return(list(df = df, 
              gbm.ATT = gbm.ATT))
}

trim_data = function(df, sample, subset_by, subset_prop, method){
  if(subset_by == "covariates"){
    if(sample == "OPEN"){
      exclude = df %>% 
        filter(S == 0, black == 0, BMI < quantile(df$BMI[which(df$S == 0 & df$black == 0)],subset_prop))
    }
    
    if(sample == "PREMIER"){
      exclude = df %>% 
        filter(S == 1, black == 1, BMI > quantile(df$BMI[which(df$S == 1 & df$black == 1)],(1-subset_prop)))
    }
  }
  
  if(subset_by == "ps"){
    data = fit_selection_model(df, method)$df
    
    if(sample == "OPEN"){
      exclude = data %>% 
        filter(S == 0, ps < quantile(data$ps[which(data$S == 0)], subset_prop))
    }
    
    if(sample == "PREMIER"){
      exclude = data %>% 
        filter(S == 1, ps > quantile(data$ps[which(data$S == 1)], (1-subset_prop)))
    }
  }
  
  df = df %>% 
    anti_join(exclude, by = "ID")
  
  return(df)
}

analysis = function(data, timepoint, method = 'gbm', trim_pctile = NULL, subset_by = NULL,sample = "OPEN", subset_prop = .05){
  if(!is.null(trim_pctile)){
    if(!is.numeric(trim_pctile)){
      stop("Trimming percentile must be numeric")
    }
    if(!(trim_pctile <= 1 & trim_pctile > 0)){
      stop("Trimming percentile must be expressed as a decimal")
    }
  }
  
  if(!is.null(subset_prop)){
    if(!is.numeric(subset_prop)){
      stop("Subsetting proportion must be numeric")
    }
    if(!(subset_prop <= 1 & subset_prop > 0)){
      stop("Subsetting proportion must be expressed as a decimal")
    }
  }
  
  ### Subset the data to the correct timepoint in the trial
  df = data %>% 
    filter(S == 0 | VISIT == timepoint) %>% 
    filter(!is.na(log_sodium_self) & !is.na(log_sodium_urine))
  
  # Subset the data according to the modified data investigation
  if(!is.null(subset_by)){
    df = trim_data(df, sample, subset_by, subset_prop, method)
  }
  
  ### Calculate mu0 in validation data unweighted
  mu0_v = with(df %>% filter(S == 0, T == 0), mean(log_sodium_self - log_sodium_urine, na.rm=TRUE))
  mu0_v = signif(mu0_v,3)
  
  ### Calculate mu0rct in the trial for the given timepoint
  mu1_rct = with(df %>% filter(S == 1, T == 1), mean(log_sodium_self - log_sodium_urine, na.rm=TRUE))
  mu0_rct = with(df %>% filter(S == 1, T == 0), mean(log_sodium_self - log_sodium_urine, na.rm=TRUE))
  
  mu0_rct = signif(mu0_rct,3)
  
  ### Fit model
  fitted_model = fit_selection_model(df, method)
  test = fitted_model$df
  test$S_lab = ifelse(test$S == 0, "Validation","Intervention Study")
  
  # Trim the weights to the specified percentile
  if(!is.null(trim_pctile)){
    test$weight[which(test$weight > quantile(test$weight[which(test$S == 0)], trim_pctile))] = quantile(test$weight[which(test$S == 0)], trim_pctile)
  }
  
  # Calculate asmd of probabilities
  ASMD = asmd(test$ps[which(test$S == 0)],test$ps[which(test$S == 1)])
  
  # Calculate weighted mu0v
  mu0_v_weighted = with(test %>% filter(S == 0, T == 0), weighted.mean((log_sodium_self - log_sodium_urine), w = weight, na.rm=TRUE))
  mu0_v_weighted = signif(mu0_v_weighted,3)
  
  # Distribution of propensity scores
  p = test %>%
    ggplot() +
    geom_density(aes(x = ps, group = as.factor(S_lab), fill = as.factor(S_lab)),alpha=.2) +
    labs(x = "Predicted Probability of Trial Membership",
         fill = "Sample") +
    ggtitle("Predicted Trial Membership Probability by Sample")
  
  # ggsave(paste0("data_example/figures/ps.png"),p,width=11,height=7,units="in",scale=1.4)
  
  # Histogram of weight values
  weight_hist = test %>% filter(S == 0) %>%
    ggplot() + geom_histogram(aes(x = weight),bins = 20) +
    ggtitle(ifelse(!is.null(trim_pctile),"Trimmed Weights in OPEN","Weights in OPEN"))
  
  # ggsave(paste0("data_example/figures/weight_hist.png"),weight_hist,width=11,height=7,units="in",scale=1.4)
  
  # Love plot
  if(method == "gbm"){
    loveplot = love.plot(cobalt::bal.tab(fitted_model$gbm.ATT))
  } else{
    loveplot = "no love plot"
  }
  
  # ggsave(paste0("data_example/figures/loveplot.png"),loveplot,width=11,height=7,units="in",scale=1.4)
  
  return(list(results = data.frame(PS_method = method,
                                   Subset_by = ifelse(is.null(subset_by),"none",subset_by),
                                   sample_subset = sample,
                                   prop_excluded = subset_prop,
                                   Method = c("Unweighted","Weighted"),
                                   mu0v = c(mu0_v,mu0_v_weighted),
                                   mu0rct = mu0_rct,
                                   bias_mu0 = c(mu0_v - mu0_rct, mu0_v_weighted - mu0_rct),
                                   ASMD = ASMD),
              prob_dist = p,
              weight_hist = weight_hist,
              loveplot = loveplot,
              data = test,
              subsetting = ifelse(is.null(subset_by), "No Modification", paste0(subset_prop*100,"% removed from ",sample))))
}
