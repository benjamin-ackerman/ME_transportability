### Setup ###
# Load Packages
library(haven);library(dplyr);library(stringr)

# Paths to OPEN and PREMIER data
premier_path = "data_example/data/PREMIER/SASv8_FILES/"
open_path = "data_example/data/OPEN/"

### Load and Clean Data ###
## OPEN 
load(paste0(open_path,"open_outliers_removed_kcal.Rdata"))

open = open %>% 
  dplyr::select(-kcal1) %>% 
  mutate(ID = as.character(id),
         S = 0, T = 0,
         log_sodium_urine = (log.urine01 + log.urine02)/2) %>% 
  dplyr::select(ID, S, T, sex = male, age_cat = agecat1, BMI = bmi, black, education, 
                log_sodium_urine = log_sodium_urine, log_sodium_self = log.self0.avg)

## PREMIER
# Note that original TX variable was: 1 - control, 2 - Tx, 3 - DASH + Tx
# Demographics
premier_dem = read_sas(paste0(premier_path,"masterp.sas7bdat"))

# Education variable
premier_ed = read_sas(paste0(premier_path,"pthistp.sas7bdat"))

# Urinary sodium (True outcome - FOR TESTING METHOD ONLY)
premier_urine = read_sas(paste0(premier_path,"labp.sas7bdat"))

# Self-reported sodium (Mismeasured outcome)
premier_selfreport = read_sas(paste0(premier_path,"foodanalp.sas7bdat"))

### Merge PREMIER self-report and sodium data
premier_sodium = premier_selfreport %>% 
  left_join(premier_urine, by = c("STUDY_ID","VISIT","TX","COHORT")) %>% 
  dplyr::select(STUDY_ID, VISIT, TX, COHORT, sodium_self = RINA, sodium_urine = MGURNA) %>% 
  mutate(log_sodium_self = log(sodium_self), log_sodium_urine = log(sodium_urine)) %>% 
  filter(VISIT == 8) %>%
  mutate(VISIT = "6mo") %>% 
  bind_rows(premier_selfreport %>%
              left_join(premier_urine, by = c("STUDY_ID","VISIT","TX","COHORT")) %>% 
              dplyr::select(STUDY_ID, VISIT, TX, COHORT, sodium_self = RINA, sodium_urine = MGURNA) %>% 
              mutate(log_sodium_self = log(sodium_self), log_sodium_urine = log(sodium_urine)) %>% 
              filter(VISIT == 10) %>% 
              mutate(VISIT = "18mo"))

### Merge all of the PREMIER data together
premier = premier_dem %>% 
  left_join(premier_ed, by = c("STUDY_ID", "COHORT","TX")) %>% 
  left_join(premier_sodium, by = c("STUDY_ID","COHORT","TX")) %>% 
  mutate(ID = str_replace(STUDY_ID, "PREMIER",""),
         S = 1, T = ifelse(TX == 1, 0, 1),
         BMI = WEIGHT_BASE/((HEIGHT/100)^2),
         sex = abs(SEX - 2)) %>% 
  dplyr::select(VISIT, ID, S, T, sex, age_cat = AGE_REL, BMI, black = RACE_REL, education = EDU_REL, 
         log_sodium_self, log_sodium_urine)

### Combined PREMIER and OPEN data, and refactor age ###
data = premier %>% 
  bind_rows(open %>% 
              mutate(VISIT = NA)) %>% 
  mutate(age_cat = forcats::fct_collapse(as.factor(age_cat),
                        `1` = c("1","2","3"),
                        `2` = "4",
                        `3` = "5",
                        `4` = "6",
                        `5` = "7",
                        `6` = c("8","9","10")))

### Save merged data ###
saveRDS(data,"data_example/data/merged_data.rds")
