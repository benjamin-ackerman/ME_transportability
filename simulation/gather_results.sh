cat calibration_sim.o* >> all.txt
Rscript get_run_times.R
rm all.txt
rm calibration_sim.o*
rm calibration_sim.e*

Rscript gather_results.R
