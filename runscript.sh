#!/bin/bash
#SBATCH -n 48             # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 7-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared           # Partition to submit to
#SBATCH --mem=64000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o myoutput_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e myerrors_%j.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH --mail-user=kjosey@hsph.harvard.edu
module load R/4.0.5-fasrc01 #Load R module
export R_LIBS_USER=$HOME/apps/R_4.0.5:$R_LIBS_USER
Rscript ~/causal-me/sim2.R