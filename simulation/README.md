# Simulation Code for "Transportability of Outcome Measurement Error Models: from Validation Studies to Intervention Trials"

This repo contains all of the necessary code to run the simulations for the manuscript titled "Transportability of Outcome Measurement Error Models: from Validation Studies to Intervention Trials." The code is set up to run on a remote cluster (originally run on the JHPCE cluster at Johns Hopkins).

### Setup

Before running the simulations on the cluster, make sure you have all of the files from this repo's `code/` directory in this `simulation/` folder, along with the following sub-directories set up:

- `pop_data/` (for the created population data files)
- `sim_results/` (for the created simulation results)
- `combined_results/` (for the combined results)

You can do this by running the following code in the command line:
```
mkdir pop_data
mkdir sim_results
mkdir combined_results
```

### config.R
The `config.R` file defines *all* of the parameters to run the simulations, including parameters you wish to vary. There are three major objects to define in the config file:

- `n_samples`: This specifies how many times you wish to draw trial/validation data from the particular binary "S" assignment. Default is set to 1.
- `pop_data_params`: This list specifies all parameters to create the population data, including
    - `scale`: how strongly the Xs impact the sample membership model
    - `b1`, `b2`, `b3` & `b4`: relative magnitudes of X coefficients on the sample membership model
    - `s_true`: list of true S model forms to simulate (each S model is specified by a list of terms to add to the main effects model)
- `sample_data_params`: This list specifies the parameters to generate the trial and validation samples and the measurement error model
    - `y_scale`: how strongly the Xs impact the measurement error model
    - `n`: sample size of the trial and sample size of the validation data
    - `g0`, `g1`, `g2` & `sigma2_z`: Parameters to simulate the true outcome Z
    - `a1`, `a2`, `a3`, `a4`, `a5`, `a6`: Parameters to simulate the mismeasured outcome Y (`a3`-`a6` are the relative magnitudes of the X coefficients)

### Running the simulations on the JHPCE cluster

1. Update and upload your `config.R` file to the cluster directory.
2. Create the population data files (this also checks if the population data already exist, in which case they're not re-generated)
```
qsub -N mkpop -cwd -l mem_free=1G,h_vmem=1G mkpop.sh 
```
3. Run the simulations (`sim.sh` submits 1000 jobs with 1000 different seeds, and each creates a results file in `sim_results/`):
```
sh sim.sh
```
4. Gather the simulation run times and combine the results across the 1000 runs (this is done in an active session. Make sure to navigate to the directory after typing `qrsh`):
```
qrsh
sh gather_results.sh
exit
```
