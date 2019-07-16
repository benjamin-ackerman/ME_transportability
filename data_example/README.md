# Data Example Code for "Transportability of Outcome Measurement Error Models: from Validation Studies to Intervention Trials"

This repo contains all of the necessary code to replicate the data example in the manuscript titled "Transportability of Outcome Measurement Error Models: from Validation Studies to Intervention Trials." 

### R scripts

- `data_cleaning.R` contains all code for pre-processing the data from the two studies. This script also merges the two studies, harmonizing and dichotomizing covariates for consistency.
- `EDA.R` contains all code for creating the exploratory analysis tables in the manuscript and Web Appendix.
- `analysis_functions.R` contains the functions used to run the analyses, fitting both the unweighted and weighted measurement error models.
- `analysis.R` executes the functions in `analysis_functions.R` and creates the table of results in the manuscript.
- `utils.R` contains miscellaneous utility functions used across the analysis stages.

### Data
Data for the validation sample, OPEN, and the lifestyle intervention trial, PREMIER, are not publicly available. 