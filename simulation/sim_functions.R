### Generate data (with multiple x's)
data_pop = function(x, scale, b1, b2, b3, b4, s_true_model){
  N = nrow(x)

  if(ncol(x) != 4){stop("x must have four columns!")}
  
  colnames(x) = paste0("x",1:4)
  
  # Get values for p5 - p9 based on true model
  s_true = s_true_params(scale, b1, b2, b3, b4, s_true_model)
  
  # return the formula for true S model
  s_true_formula = s_formula(s_true_model)
  
  # Specify the true propensity score model
  p = expit(scale*(b1*x[,1] + b2*x[,2] + b3*x[,3] + b4*x[,4] + 
                     s_true$b5*x[,3]^2 + s_true$b6*x[,4]^2 + s_true$b7*x[,3]*x[,4] + s_true$b8*x[,1]*x[,4] + s_true$b9*x[,1]*x[,3]))
  
  # Gather the data together
  data = data.frame(p = p, 
                    x1 = x[,1], x2=x[,2], x3=x[,3],x4=x[,4],
                    stringsAsFactors = FALSE)
  
  return(list(data = data,
              scale = scale, 
              b1 = b1, b2 = b2, b3 = b3, b4 = b4, 
              b5 = s_true$b5, b6 = s_true$b6, b7 = s_true$b7, b8 = s_true$b8, b9 = s_true$b9,
              s_true = s_true_formula))
}

assign_S = function(df_pop_obj){
  df_pop = df_pop_obj$data
  
  N = nrow(df_pop)
  
  # Generate S from the sample selection probabilities and then randomly sample the trial and validation data
  df_pop$S = rbinom(N, 1, prob = df_pop$p)
  df_pop$S_lab = factor(df_pop$S, levels=c(0,1),labels=c("Validation","Trial"))
  
  df_pop_obj$data = df_pop
  
  return(df_pop_obj)
}
  
data_samples = function(df_pop_obj, n, g0, g1, g2, sigma2_z, y_scale, a1, a2, a3, a4, a5, a6,seed = NULL){
  if(!is.null(seed)){set.seed(seed)}

  df_pop = df_pop_obj$data

  n_v = n_rct = n
  
  # Gather the data together
  data = df_pop %>%
    group_by(S) %>%
    do(sample_n(.,n_rct)) %>% 
    ungroup()
  
  # Treatment assignment
  data$T = c(rep(0, n),rbinom(n, 1, .5))
  
  # Simulate Z's
  data$Z0 = rnorm(n_v + n_rct,g0 + g2,sqrt(sigma2_z))
  data$Z1 = rnorm(n_v + n_rct,g0 + g1 + g2,sqrt(sigma2_z))
  
  # Get rid of random noise in Z (there's an issue with X here though! They can't be collinear)
  #Z0 = g0 + g2*X
  #Z1 = g0 + g1 + (g2+g3)*X
  
  data$Z = (1-data$T)*data$Z0 + data$T*data$Z1
  
  #Y = a_1*Z + a_2*T + a_3*X1 + a_4*X2 + a_5*X3 + a_6*X4 
  
  # Scenarios we care about: 
  # all four impact S but not Y (a3:a6 = 0, p1:p4 =/= 0)
  # all four impact Y but not S (a3:a6 =/= 0, p1:p4 = 0)
  # two impact S, two impact Y (a3,a4,p3,p4 = 0, a5,a6,p1,p2 =/= 0)
  # two impact both, two only impact S
  # two impact both, two only impact Y
  
  # Center x's so that they're mean 0 in the trial:
  #data$x1 = data$x1 - mean(data$x1[which(data$S == 1)])
  #data$x2 = data$x2 - mean(data$x2[which(data$S == 1)])
  #data$x3 = data$x3 - mean(data$x3[which(data$S == 1)])
  #data$x4 = data$x4 - mean(data$x4[which(data$S == 1)])
  
  data$Y0 = rnorm(n_v + n_rct, a1*data$Z0 + y_scale*(a3*data$x1 + a4*data$x2 + a5*data$x3 + a6*data$x4), sqrt(sigma2_z*1.5))
  data$Y1 = rnorm(n_v + n_rct, a1*data$Z1 + a2 + y_scale*(a3*data$x1 + a4*data$x2 + a5*data$x3 + a6*data$x4), sqrt(sigma2_z*1.5))
  
  data$Y = (1-data$T)*data$Y0 + data$T*data$Y1
  
  return(list(data = data,
              covariate_diffs = x_means(data),
              scale = df_pop_obj$scale, 
              b1 = df_pop_obj$b1, b2 = df_pop_obj$b2, b3 = df_pop_obj$b3, b4 = df_pop_obj$b4, 
              b5 = df_pop_obj$b5, b6 = df_pop_obj$b6, b7 = df_pop_obj$b7, b8 = df_pop_obj$b8, b9 = df_pop_obj$b9,
              s_true = df_pop_obj$s_true,
              y_scale = y_scale,
              a3 = a3, a4 = a4, a5 = a5, a6 = a6))
}

data_gen = function(pop, #population data created from data_pop (saved as RDS)
                    n_samples, # Number of times to randomly sample from each population
                    n, # a list of quadratic/interaction terms to include in model (blank list if simple)
                    g0, g1, g2, sigma2_z, # Z ~ T + X model
                    y_scale, a1, a2, a3, a4, a5, a6){
  
  pop_S = assign_S(pop)

  if(n_samples == 1){
    samples = data_samples(pop_S, n, g0, g1, g2, sigma2_z, y_scale, a1, a2, a3, a4, a5, a6, NULL)
  } else{
    samples = 1:n_samples %>% 
      map(~data_samples(pop_S, n, g0, g1, g2, sigma2_z, y_scale, a1, a2, a3, a4, a5, a6, .))
  }
  
  return(samples)
}

### Specify the weighting model
weight_model = function(model, data_obj){
  # model options are "unweighted", "weighted_x", "weighted_x2", "weighted_y" or "weighted_y2"
  dat = data_obj$data
  
  if(model == "unweighted"){
    ps = 0.5
  }
  if(model == "simple"){
    ps = predict(glm(as.formula(s_formula(list())), data = dat,family='binomial'),type='response')
  }
  if(model == "true"){
    ps = predict(glm(as.formula(data_obj$s_true), data = dat,family='binomial'),type='response')
  }
  return(ps)
}

### Create scenarios for true selection model
s_true_params = function(scale, b_1, b_2, b_3, b_4, terms){
  b5 = 0; b6 = 0; b7 = 0; b8 = 0; b9 = 0
  if("x3^2" %in% unlist(terms)){
    b5 = 0.5*b_3
  }
  if("x4^2" %in% unlist(terms)){
    b6 = 0.5*b_4
  }
  if("x3*x4" %in% unlist(terms)){
   b7 = 0.5*(b_3 + b_4)/2
  }
  if("x1*x4" %in% unlist(terms)){
    b8 = 0.5*(b_1 + b_4)/2
  }
  if("x1*x3" %in% unlist(terms)){
    b9 = 0.5*(b_1 + b_3)/2
  }
  return(list(b5 = b5,b6=b6,b7=b7,b8=b8,b9=b9))
}

### Turn list of parameters to include in model to a formula
s_formula = function(s_true_model){
 model = "S ~ -1 + x1 + x2 + x3 + x4" 
 
 vars = unlist(s_true_model)
 # Get rid of any fake variables that shouldn't belong
 vars = vars[stringr::str_detect(vars,"^x[1-4]")]
 vars = stringr::str_replace(vars,"_","^")
 vars = stringr::str_replace(vars,"-","*")
 
 if(length(vars) == 0){
   return(model)
 }
 else{
   model = paste0(model," + ", paste("I(",vars,")",sep="",collapse=" + "))
   return(model)
 }
}

# Look at the mean difference between the validation and trial
x_means = function(data){
  means = data %>%
    group_by(S) %>% 
    summarise_at(.vars=c("x1","x2","x3","x4"), .funs=mean) %>% 
    dplyr::select(-S) %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(diff = V1-V2,
           covariate = row.names(.)) %>% 
    dplyr::select(covariate,diff)
  
  return(means)
}

# Calculate Error under control in data set s using simulated data
get_mu = function(data_obj, s, weighting){
  # Extract data frame from data generating object
  dat = data_obj$data
  
  ### Get the predicted sample membership probabilities and generate weights
  ps = weight_model(weighting, data_obj)
  
  dat$weights = ifelse(dat$S == 0, ps/(1-ps), 1)
  
  mu = with(dat %>% filter(S == s, T == 0), weighted.mean(Y - Z,weights))
  
  return(mu)
}

# Gather and return all results
get_results = function(data_obj, weighting){
  ### DoM:
  # Propensity scores under true and specified models
  prob_true = weight_model("true",data_obj)
  prob_spec = weight_model("simple", data_obj)
  
  # Calculate the asmd of the true propensity scores
  ASMD = asmd(prob_true[which(data_obj$data$S == 0)], prob_true[which(data_obj$data$S == 1)])
  
  # Calculate the degree of misspecification
  dom = DoM(prob_true, prob_spec)
  
  # Calculate mu_v
  mu_v = get_mu(data_obj, 0, weighting)

  # Calculate mu_l
  # Try fitting the model ONLY in the validation data *ALWAYS* with the specified model
  mu_l = get_mu(data_obj, 1, "unweighted")

  bias_mu0 = mu_v - mu_l
  
  results = data.frame(s_true = data_obj$s_true, 
                       weighting = weighting,
                       scale = data_obj$scale,
                       y_scale = data_obj$y_scale,
                       DoM = dom,
                       ASMD = ASMD,
                       m_01_bias = bias_mu0, 
                       mu_v = mu_v,
                       mu_l = mu_l,
                       stringsAsFactors = FALSE)
  return(results)
}

### Run the simulations using the population data by specifying the seed (to sample 1000 times...)
run_samples = function(df_pop_obj, dat_params, seed){
  data_objs = seq_along(dat_params) %>% 
    map(~data_samples(
      df_pop_obj, dat_params[[.]]$n, dat_params[[.]]$g0, dat_params[[.]]$g1, 
      dat_params[[.]]$g2, dat_params[[.]]$sigma2_z, dat_params[[.]]$y_scale, 
      dat_params[[.]]$a1, dat_params[[.]]$a2, dat_params[[.]]$a3, dat_params[[.]]$a4, 
      dat_params[[.]]$a5, dat_params[[.]]$a6, seed
    ))
  
  fit_params = list(weighting = c("unweighted","simple","true"),
                    data_obj = data_objs) %>%
    cross()
  
  results = seq_along(fit_params) %>% 
    map_df(~get_results(fit_params[[.]]$data_obj, fit_params[[.]]$weighting))
  return(results)
}

run_pops = function(population, dat_params, n_samples){
  pop = assign_S(readRDS(paste0("pop_data/pop_",population,".rds")))
  
  if(n_samples == 1){
    results = run_samples(pop, dat_params, NULL)
  } else{
    results = 1:n_samples %>% 
      map_df(~run_samples(pop, dat_params, .)) %>% 
      group_by(s_true,weighting, scale,y_scale) %>% 
      summarise(m_01_bias = mean(m_01_bias),
                mu_v = mean(mu_v),
                mu_l = mean(mu_l),
                DoM = mean(DoM),
                ASMD = mean(ASMD)) %>% 
      ungroup()
  }
  
  return(results)
}

##### PLOTTING FUNCTIONS #####
plot_bias = function(data, s_model = NULL){
  p_dat = data %>% 
    filter(weighting == "simple" & scale > 0 & y_scale == "scale[y] == 0.2") %>% 
    mutate(s_true = forcats::fct_reorder(s_true,DoM))
  
  p = data %>% ungroup() %>% 
    mutate(s_true = forcats::fct_relevel(s_true, levels(p_dat$s_true))) %>% 
    ggplot() +
    geom_line(aes(x = ASMD, y = abs(bias), col=weighting)) +
    facet_grid(y_scale ~ s_true, label = 'label_parsed') +
    labs(color = "Weighting Method",
         x = "ASMD of true selection probabilities",
         y = expression(paste({'|Bias'}[mu[0]],'|'))) +
    ggtitle(expression(paste('Absolute Bias of estimating ',mu[0], ' in the trial using validation data')))
  return(p)
}

plot_bias_by_dom = function(data){
  d1 = data %>% 
    filter(weighting == "simple") %>% 
    dplyr::select(s_true, scale, y_scale, DoM)
  
  d2 = data %>% 
    dplyr::select(-DoM) %>% 
    left_join(d1, by = c("s_true","scale","y_scale")) %>% 
    mutate(scale = paste0("scale[S] == ",scale))
  
  p_dat = d2 %>% 
    filter(weighting == "simple" & scale > 0 & y_scale == "scale[y] == 0.2") %>% 
    mutate(s_true = forcats::fct_reorder(s_true,DoM))
  
  p = d2 %>% ungroup() %>% 
    mutate(s_true = forcats::fct_relevel(s_true, levels(p_dat$s_true))) %>% 
    ggplot() +
    geom_line(aes(x = DoM, y = abs(bias), col=weighting)) +
    facet_grid(y_scale ~ scale, label = 'label_parsed') +
    labs(color = "Weighting Method",
         x = "DoM",
         y = expression(paste({'|Bias'}[mu[0]],'|'))) +
    ggtitle(expression(paste('Absolute Bias of estimating ',mu[0], ' in the trial using validation data')))
  return(p)
}

plot_scale_to_asmd = function(data, s_model = NULL){
  p = ggplot(results) + geom_line(aes(x = scale,y = ASMD, group = s_true, color = s_true)) +
    scale_colour_discrete(labels = scales::parse_format()) +
    labs(color = "Specified Selection Model", title = "True selection model scale parameter vs. ASMD")
  return(p)
}

plot_dom = function(data){
  p_dat = data %>% 
    filter(weighting == "simple" & scale > 0 & y_scale == "scale[y] == 0.2") %>% 
    mutate(s_true = forcats::fct_reorder(s_true,DoM))
  
  p = p_dat %>% ggplot() + 
    geom_line(aes(x = s_true,y = DoM, group = as.factor(scale), col=as.factor(scale))) +
    scale_x_discrete(labels = new_parse_format(levels(p_dat$s_true)))+
    labs(x = "True S model form",
         col = "Scale") +
    theme(axis.text.x = element_text(size=6)) + 
    ggtitle("Degree of Misspecification by True S Model Form and Parameter Scale")
  
  return(p)
}


