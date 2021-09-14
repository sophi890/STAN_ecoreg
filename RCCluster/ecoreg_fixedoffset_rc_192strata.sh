#!/bin/bash

#SBATCH -n 48 #number of cores
#SBATCH -t 6-23:59
#SBATCH -o /n/home07/swoodward/STAN_ecoreg/test3_error.out #specify where to save errors returned by the program
#SBATCH -e /n/home07/swoodward/STAN_ecoreg/test3_log.err #specify where to save the output log
#SBATCH --array=1 #number of jobs to run, it is currently set to 1 job(1 year), change it to array=1-13 for 13 years(jobs)
#SBATCH --mem=80000 #memory requested
#SBATCH -J fit192.fixedoffsets  #job name
#SBATCH --mail-type=END #notifications for job done
#SBATCH --mail-user=swoodward@college.harvard.edu # send to address

export R_LIBS_USER=/n/home07/swoodward/apps/R_3.5.1:$R_LIBS_USER  #change THIS accordingly
module load gcc/9.2.0-fasrc01 R/3.6.3-fasrc02
#module load R/3.6.2-fasrc01
#module load R_core/3.4.2-fasrc01
#module load R_packages/3.4.2-fasrc02

R CMD BATCH --quiet --no-restore --no-save ecoreg_fixedoffset_rc_192strata.R /n/home07/swoodward/STAN_ecoreg/test3_${SLURM_ARRAY_TASK_ID}.Rout
