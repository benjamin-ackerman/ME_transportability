### Expit function for propensity scores
expit = function(x){
  return(exp(x)/(1+exp(x)))
}

### Calculate the Degree of Misspecification
DoM = function(pred_true, pred_spec){
  sigma2_true = sqrt(var(pred_true))
  dom = sum(abs(pred_spec - pred_true))/(sigma2_true*length(pred_true))
  
  return(dom)
}

### Calculate the ASMD
asmd = function(x1, x2){
  pooled_sd = sqrt(((length(x1)-1)*var(x1) + (length(x2)-1)*var(x2))/(length(x1)+length(x2)-2))
  return(abs((mean(x1)-mean(x2))/pooled_sd))
}

### Calculate coverage
get_weighted_se = function(Y, Z, weights){
  n = length(Y)
  
  X = cbind(Y,Z)
  W = as.vector(weights)
  
  vcov = cov.wt(X, W)$cov
  
  vpooled = vcov[1,1] + vcov[2,2] - 2*vcov[1,2]
  se = sqrt(vpooled/n)
  
  return(se)
}

coverage = function(estimate, se, truth){
  CI_lower = estimate - 1.96*se
  CI_upper = estimate + 1.96*se
  
  covered = ifelse(truth >= CI_lower & truth <= CI_upper, 1, 0)
  
  return(covered)
}

### Function from github help to parse axis labels as equations:
new_parse_format <- function(text) {
  text <- as.character(text)
  out <- vector("expression", length(text))
  for (i in seq_along(text)) {
    expr <- parse(text = text[[i]])
    out[[i]] <- if (length(expr) == 0)
      NA
    else expr[[1]]
  }
  out
}
