#!/bin/bash
#SBATCH --job-name pov-sims
#SBATCH --partition short
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 3
#SBATCH --time 0-4:00
#SBATCH --mem-per-cpu=500
#SBATCH --array=1-50
#SBATCH -o myscript_%j.out
#SBATCH -e myscript_%j.err
#SBATCH --mail-type=NONE              # Type of email notification- NONE,BEGIN,END,FAIL,ALL
#SBATCH --mail-user=rritaz@uw.edu
source /etc/profile.d/z00_lmod.sh
module load R
Rscript fund-alloc.R $SLURM_ARRAY_TASK_ID
