### Setup
library(dplyr);library(generalize)

### Load data
data = readRDS("data_example/data/merged_data.rds")
data$age_cat = as.factor(data$age_cat)
data$education = as.factor(data$education)

n = data %>% group_by(S) %>%filter(!duplicated(ID)) %>%  summarise(n()) %>% pull()

### Covariate table
covtab = covariate_table("S",names(data)[5:9], data %>% group_by(S) %>% filter(!duplicated(ID)) %>% ungroup(),weighted_table = FALSE)
colnames(covtab) = c(paste0("OPEN (n=",n[1],")"),paste0("PREMIER (n=",n[2],")"),"ASMD")
rownames(covtab) = c("Male", "<= 40", "41-45","46-50","51-55","56-60",">= 61","BMI","Black","College","Grad School")

### Table of outcome means by treatment
outcomes = data %>% 
  group_by(S,VISIT, T) %>% 
  summarise(self = mean(log_sodium_self,na.rm=TRUE),
            urine = mean(log_sodium_urine, na.rm=TRUE),
            epsilon = self - urine) %>% ungroup() %>%
  filter(!(S == 1 & is.na(VISIT))) %>% 
  dplyr::select(-S,-T,-VISIT) %>% 
  t() %>% as.data.frame()

outcomes = outcomes[,c(1,4,5,2,3)]
names(outcomes) = c("Control","Control","Treatment","Control","Treatment")
row.names(outcomes) = c("Self-Reported (Y)","Urine (Z)","$\\mu_{a}^{s}$")

error_tx = function(data, timepoint){
  df = data %>% 
    filter(S == 0 | VISIT == timepoint) %>% 
    filter(!is.na(log_sodium_self) & !is.na(log_sodium_urine))
  
  error_treatment = with(df %>% filter(S == 1), t.test((log_sodium_self - log_sodium_urine) ~ T))
  error_treatment = broom::tidy(error_treatment)
  error_treatment$'95% CI' = paste0("(",signif(error_treatment$conf.low,3),", ",signif(error_treatment$conf.high,3),")")
  error_treatment = error_treatment[c("estimate","95% CI","p.value")]
  names(error_treatment) = c("$\\text{bias}_{ATE}$","$95\\%$ CI","p-value")
  
  return(error_treatment)
}

error_treatment = cbind(Timepoint = c("18 months","6 months"),
                        error_tx(data,"18mo") %>%
                          bind_rows(error_tx(data,"6mo")))
error_treatment = error_treatment[c(2,1),]
row.names(error_treatment) = NULL

# Fit a measurement error model to the full data using the covariates
model_covariates = lm(log_sodium_self ~ log_sodium_urine + sex+ age_cat + BMI + as.factor(education) + black, 
   data = data %>% filter(T == 0,S == 0))
summary(model_covariates)
