#######################################################################
# outdated functions #
#######################################################################
data_gen = function(n, x, # true Y model, sample sizes
                    scale, b1, b2, b3, b4, 
                    s_true_model, # a list of quadratic/interaction terms to include in model (blank list if simple)
                    g0, g1, g2, sigma2_z, # Z ~ T + X model
                    a1, a2, a3, a4, a5, a6){ # Y ~ Z + X model
  N = nrow(x)
  n_v = n_rct = n
  
  if(ncol(x) != 4){stop("x must have four columns!")}
  
  colnames(x) = paste0("x",1:4)
  
  # Get values for p5 - p9 based on true model
  s_true = s_true_params(b1, b2, b3, b4, s_true_model)
  
  # return the formula for true S model
  s_true_formula = s_formula(s_true_model)
  
  # Specify the true propensity score model
  p = expit(scale*(b1*x[,1] + b2*x[,2] + b3*x[,3] + b4*x[,4] + 
                     s_true$b5*x[,3]^2 + s_true$b6*x[,4]^2 + s_true$b7*x[,3]*x[,4] + s_true$b8*x[,1]*x[,4] + s_true$b9*x[,1]*x[,3]) + rnorm(N))
  
  # Generate S from the sample selection probabilities and then randomly sample the trial and validation data
  S = rbinom(N, 1, prob = p)
  
  # Gather the data together
  data = data.frame(S = S,
                    p = p, 
                    x1 = x[,1], x2=x[,2], x3=x[,3],x4=x[,4],
                    S_lab = factor(S, levels=c(0,1),labels=c("Validation","Trial")), 
                    stringsAsFactors = FALSE) %>% 
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
  
  data$Y0 = rnorm(n_v + n_rct, a1*data$Z0 + a3*data$x1 + a4*data$x2 + a5*data$x3 + a6*data$x4, sqrt(sigma2_z*1.5))
  data$Y1 = rnorm(n_v + n_rct, a1*data$Z1 + a2 + a3*data$x1 + a4*data$x2 + a5*data$x3 + a6*data$x4, sqrt(sigma2_z*1.5))
  
  data$Y = (1-data$T)*data$Y0 + data$T*data$Y1
  
  return(list(data = data,
              covariate_diffs = x_means(data),
              scale = scale, 
              b1 = b1, b2 = b2, b3 = b3, b4 = b4, 
              b5 = s_true$b5, b6 = s_true$b6, b7 = s_true$b7, b8 = s_true$b8, b9 = s_true$b9,
              s_true = s_true_formula,
              a3 = a3, a4 = a4, a5 = a5, a6 = a6))
}

### OLD: Generate a data set
data_gen = function(y_true, n_trial, n_calib, # true Y model, sample sizes
                     cov, b0, b1, # X
                     g0, g1, g2, g3, sigma2_z, # Z
                     a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10){ # Y
   N = n_trial + n_calib
   S = c(rep(1,n_trial),rep(0,n_calib))

   #generate two rvs with mean 0 and sd 1 that are correlated
   XU = mvrnorm(N, mu = c(0,0), Sigma = matrix(c(1,cov,cov,1),nrow=2))

   #Grab the first column of XU to be X, and add a0 + a1*S
   X = XU[,1] + b0 + b1*(S==0)
   U = XU[,2] #For now, don't do anything with U, but can specify later

   # Treatment assignment
   T = c(rbinom(n_trial, 1, .5), rep(0, n_calib))

   # Simulate Z's
   Z0 = rnorm(N,g0 + g2*X,sqrt(sigma2_z))
   Z1 = rnorm(N,g0 + g1 + (g2+g3)*X,sqrt(sigma2_z))

   # Get rid of random noise in Z (there's an issue with X here though! They can't be collinear)
   #Z0 = g0 + g2*X
   #Z1 = g0 + g1 + (g2+g3)*X

   Z = (1-T)*Z0 + T*Z1
   # Simulate Y's according to true underlying measurement error model
   if(y_true == "x"){
     Y0 = rnorm(N, a0 + a1*Z0 + a3*X, sqrt(sigma2_z*1.5))
     Y1 = rnorm(N, a0 + a1*Z1 + a2 + a3*X, sqrt(sigma2_z*1.5))

     a4 = a5 = a6 = a7 = a8 = a9 = a10 = 0
   }

   if(y_true == "x2"){
     Y0 = rnorm(N, a0 +a1*Z0 + a3*X + a5*(X^2), sqrt(sigma2_z*1.5))
     Y1 = rnorm(N, a0 +a1*Z1 + a2 + a3*X + a5*(X^2), sqrt(sigma2_z*1.5))

     a4 = a6 = a7 = a8 = a9 = a10 = 0
   }

   if(y_true == "z-x"){
     Y0 = rnorm(N, a0 +a1*Z0 + a3*X + a4*X*Z0, sqrt(sigma2_z*1.5))
     Y1 = rnorm(N, a0 +a1*Z1 + a2 + (a3+a7)*X + (a4+a8)*X*Z1, sqrt(sigma2_z*1.5))

     a5 = a6 = a9 = a10 = 0
   }

   if(y_true == "z-x2"){
     Y0 = rnorm(N, a0 +a1*Z0 + a3*X + a4*X*Z0 + a5*(X^2), sqrt(sigma2_z*1.5))
     Y1 = rnorm(N, a0 +a1*Z1 + a2 + (a3+a7)*X + (a4+a8)*X*Z1 + (a5+a9)*(X^2), sqrt(sigma2_z*1.5))

     a6 = a10 = 0
   }

   if(y_true == "z-x2-full"){
     Y0 = rnorm(N, a0 + a1*Z0 + a3*X + a4*X*Z0 + a5*(X^2) + a6*Z0*(X^2), sqrt(sigma2_z*1.5))
     Y1 = rnorm(N, a0 + a1*Z1 + a2 + (a3+a7)*X + (a4+a8)*X*Z1 + (a5+a9)*(X^2) +
                  (a6+a10)*Z1*(X^2), sqrt(sigma2_z*1.5))
   }

   Y = (1-T)*Y0 + T*Y1

   # Gather the data
   data = data.frame(S,X,T,Z0,Z1,Z,Y0,Y1,Y)

   # Calculate mu_01 under true Y model
   m_01_truth = a0 + a3*b0 + a5*(b0^2) + (g0 + g2*b0)*(a1 + a4*b0 + a6*(b0^2) - 1)

   return(list(data = data,
               y_true = y_true,
               m_01_truth = m_01_truth,
               b0 = b0,
               b1 = b1,
               a3 = a3,
               g0 = g0,
               g2 = g2))
 }
 ### Change sigma2_y to be 1.5 sigma2_z
# Calculate the bias
bias_calc = function(data_obj, y_spec, weighting){
   # Extract data frame from data generating object
   data = data_obj$data
#   ### Get the predicted sample membership probabilities and generate weights
   ps = weight_model(weighting, data)
   data$weights = ifelse(data$S == 0, ps/(1-ps), 1)

   ## Predictions under correct and specified models:
   pred_true = lm(y_formula(data_obj$y_true), data = data %>% filter(S == 0, T == 0))$fitted
   pred_spec = lm(y_formula(y_spec), data = data %>% filter(S == 0, T == 0))$fitted

   # Calculate the degree of misspecification
   dom = DoM(pred_true, pred_spec)

   ## Fit specified model with selected weighting option
   spec_model = lm(y_formula(y_spec), weights = weights, data = data %>% filter(S == 0, T == 0))

   ## Gather coefficients from (mis)specified model
   a = summary(spec_model)$coefficients[,1]

   a0_hat = as.numeric(a["(Intercept)"])
   a1_hat = as.numeric(a["Z"])
   a3_hat = as.numeric(ifelse(is.na(a["X"]), 0, a["X"]))
   a4_hat = as.numeric(ifelse(is.na(a["Z:X"]), 0, a["Z:X"]))
   a5_hat = as.numeric(ifelse(is.na(a["I(X^2)"]), 0, a["I(X^2)"]))
   a6_hat = as.numeric(ifelse(is.na(a["Z:I(X^2)"]), 0, a["Z:I(X^2)"]))

   ## Calculate means of Z and X (weighted if we're using weighting methods)
   # UPDATE: I don't think I need the if statements, weighted.mean will handle weights of 1
   # if(weighting == "unweighted"){
   #   ex = with(data %>% filter(S == 0), mean(X))
   #   ez = with(data %>% filter(S == 0, T == 0),mean(Z))
   # }
   #
   #  if(weighting != "unweighted"){
   ex_0 = with(data %>% filter(S == 0), weighted.mean(X,weights))
   ex_1 = with(data %>% filter(S == 1), mean(X))
   ez_0 = with(data %>% filter(S == 0, T == 0),weighted.mean(Z,weights))
   ez_1 = with(data %>% filter(S == 1, T == 0), mean(Z))
   ezx_0 = with(data %>% filter(S == 0, T == 0),weighted.mean(X*Z,weights))
   ezx_1 = with(data %>% filter(S == 1, T == 0), mean(X*Z))
   ex2_0 = with(data %>% filter(S == 0, T == 0),weighted.mean(X^2,weights))
   ex2_1 = with(data %>% filter(S == 1, T == 0), mean(X^2))
   ezx2_0 = with(data %>% filter(S == 0, T == 0),weighted.mean(Z*X^2,weights))
   ezx2_1 = with(data %>% filter(S == 1, T == 0), mean(Z*X^2))
   #}

   ##### mu_00 Estimate #####
   # Estimate mu_00 under specified Y model
   #m_00_estimate = a0_hat + a3_hat*ex + a5_hat*(ex^2) + ez*(a1_hat + a4_hat*ex +
   #                                                           a6_hat*(ex^2) - 1)

   # Find bias between mu_00 estimate and mu_01 truth
   #m_01_bias = m_00_estimate - data_obj$m_01_truth
   m_01_bias = (a1_hat - 1)*(ez_0-ez_1)+
     a3_hat*(ex_0 - ex_1) +
     a4_hat*(ezx_0 - ezx_1) +
     a5_hat*(ex2_0 - ex2_1) +
     a6_hat*(ezx2_0 - ezx2_1)

   results = data.frame(y_true = data_obj$y_true,
                        y_spec = y_spec,
                        b0 = data_obj$b0,
                        g0 = data_obj$g0,
                        b1 = data_obj$b1,
                        a3 = data_obj$a3,
                        weighting = weighting,
                        DoM = dom,
                        m_01_bias = m_01_bias,
                        #m_00_estimate = m_00_estimate,
                        #m_01_truth = data_obj$m_01_truth,
                        stringsAsFactors = FALSE)

   return(results)
 }

### Create the formula based on the model
y_formula = function(spec){
   if(spec == "simple"){
     model = as.formula("Y ~ Z")
   }

   if(spec == "x"){
     model = as.formula("Y ~ Z + X")
   }

   if(spec == "x2"){
     model = as.formula("Y ~ Z + X + I(X^2)")
   }

   if(spec == "z-x"){
     model = as.formula("Y ~ Z*X ")
   }

   if(spec == "z-x2"){
     model = as.formula("Y ~ Z*X + I(X^2)")
   }

   if(spec == "z-x2-full"){
     model = as.formula("Y ~ Z*X + Z*I(X^2)")
   }

   return(model)
 }

# Old mu function
get_mu = function(data_obj, s, y_spec, weighting){
   # Extract data frame from data generating object
   dat = data_obj$data

   ### Get the predicted sample membership probabilities and generate weights
   ps = weight_model(weighting, dat)

   dat$weights = ifelse(dat$S == 0, ps/(1-ps), 1)

   fit = lm(y_formula(y_spec), weights = weights, data = dat %>% filter(S == s, T == 0))

   a = summary(fit)$coefficients[,1]

   a0_hat = as.numeric(a["(Intercept)"])
   a1_hat = as.numeric(a["Z"])
   a3_hat = as.numeric(ifelse(is.na(a["X"]), 0, a["X"]))
   a4_hat = as.numeric(ifelse(is.na(a["Z:X"]), 0, a["Z:X"]))
   a5_hat = as.numeric(ifelse(is.na(a["I(X^2)"]), 0, a["I(X^2)"]))
   a6_hat = as.numeric(ifelse(is.na(a["Z:I(X^2)"]), 0, a["Z:I(X^2)"]))

   ex = with(dat %>% filter(S == s, T == 0), weighted.mean(X, weights))
   ez = with(dat %>% filter(S == s, T == 0),weighted.mean(Z, weights))
   ex2 = with(dat %>% filter(S ==s, T == 0), weighted.mean(X^2, weights))
   ezx = with(dat %>% filter(S == s, T == 0), weighted.mean(X*Z, weights))
   ezx2 = with(dat %>% filter(S == s, T == 0), weighted.mean(Z*X^2, weights))

   mu = a0_hat + (a1_hat-1)*ez + a3_hat*ex +a4_hat*ezx + a5_hat*ex2 + a6_hat*ezx2

   return(mu)
 }

### Plot the bias results by model specification and b0 and g0 values
plot_bias = function(data, spec = NULL, b_0 = 0, g_0 = 0, g_2 = 0){
   if(is.null(spec)){spec = unique(results$y_spec)}
 
   p = data %>% 
     filter(y_spec %in% spec, b0 == b_0, g0 == g_0, g2 == g_2) %>% #, a3 %in% c(0,0.5, 1)) %>% 
     ggplot() + 
     geom_line(aes(x = b1, y = abs(bias), group = interaction(y_spec_new, weighting), linetype = weighting, col=y_spec_new)) + 
     facet_grid(a3_new ~ y_true_new, label = "label_parsed") +
     labs(x = expression(paste(beta[1],': X diff. b/w Trial and Calibration')),
          y = expression(paste({'|Bias'}[mu[0]],"|")),
          col = "Specified Y model") +
     scale_colour_discrete(labels = scales::parse_format())+
     ggtitle(expression(paste('Absolute Bias of estimating ',mu[0], ' in the trial using validation data')))
   
   return(p)
 }

### Plot the DoM by model specification and b0 and g0 values
plot_dom = function(data, spec = NULL, b_0 = 0, g_0 = 0, g_2 = 0){
  if(is.null(spec)){spec = unique(results$y_spec)}
  
  p =  ggplot(data %>% filter(y_spec %in% spec, b0 == b_0, g0 == g_0, g2 == g_2)) + 
    geom_line(aes(x = b1, y = DoM, group = interaction(y_spec_new, weighting), linetype = weighting, col=y_spec_new)) + 
    facet_grid(a3_new ~ y_true_new, label = "label_parsed") +
    labs(x = expression(paste(beta[1],': X diff. b/w Trial and Calibration')),
         y = 'Degree of Model Misspecification (DoM)',
         col = "Specified Y model") +
    scale_colour_discrete(labels = scales::parse_format())+
    ggtitle('DoM by true Y model')
  
  return(p)
}

run = function(seed){
  set.seed(seed)
  data_objs = seq_along(pops) %>% 
    map(~data_samples(pops[[.]], 1000, 0, 2, 0, 1,
                      1, 0, 0, 1, 2, 0.5))
  
  fit_params = list(weighting = c("unweighted","simple","true"),
                    data_obj = data_objs) %>%
    cross()
  
  results = seq_along(fit_params) %>% 
    map_df(~get_results(fit_params[[.]]$data_obj, fit_params[[.]]$weighting))
  
  return(results)
}

### Delete unnecessary simulation runs
detect_extras = function(x){
  condition_1 = fit_params[[x]]$data_obj$y_true == "x" & fit_params[[x]]$y_spec %in% c("x2","z-x","z-x2","z-x2-full")
  condition_2 = fit_params[[x]]$data_obj$y_true == "x2" & fit_params[[x]]$y_spec %in% c("z-x2","z-x2-full")
  condition_3 = fit_params[[x]]$data_obj$y_true == "z-x" & fit_params[[x]]$y_spec %in% c("z-x2","z-x2-full")
  condition_4 = fit_params[[x]]$data_obj$y_true == "z-x2" & fit_params[[x]]$y_spec %in% c("z-x2-full")
  
  return(ifelse(sum(condition_1,condition_2,condition_3,condition_4) > 0,x,NA))
}
