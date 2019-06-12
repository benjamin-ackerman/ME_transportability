for f in {10001..11000}
do
qsub -N calibration_sim -cwd -l mem_free=1G,h_vmem=1G run.sh  $f
done
