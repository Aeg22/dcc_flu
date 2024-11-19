#!/bin/bash  
#SBATCH -J cellranger_aggr_degregori16
#SBATCH -o cellranger_aggr_degregori16_%j.out  # Output file with the job ID
#SBATCH -e cellranger_aggr_degregori16_%j.err  # Error file with the job ID
#SBATCH -t 1-00:00:00   # Set the wall time: D-HH:MM:SS
#SBATCH -n 1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G

# Paths
CSV=aggregation.csv

toolpath=/projects/mboorgula@xsede.org/software/cellranger-7.1.0
$toolpath/cellranger aggr \
--id=Aggr_noNorm_DeGregori_16samples \
--csv=${CSV} \
--normalize=none

