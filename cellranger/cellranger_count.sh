#!/bin/bash  
#SBATCH -J cellranger_degregori16
#SBATCH -o cellranger_degregori16_%j.out  # Output file with the job ID
#SBATCH -e cellranger_degregori16_%j.err  # Error file with the job ID
#SBATCH -t 1-00:00:00   # Set the wall time: D-HH:MM:SS
#SBATCH -n 1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G

toolpath=/projects/mboorgula@xsede.org/software/cellranger-7.1.0
$toolpath/cellranger multi --id=output_multi --csv=file.csv
