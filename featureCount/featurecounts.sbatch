#!/bin/bash
#SBATCH --job-name=<job_name>				 # Job name
#SBATCH --mail-type=FAIL,END				 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<email>				 # Where to send mail
#SBATCH --nodes=1					 # Number of cores job will run on
#SBATCH --ntasks=1					 # Number of CPU (processers, tasks)
#SBATCH --time=01:00:00					 # Time limit hrs:min:sec
#SBATCH --mem=4gb					 # Memory limit
#SBATCH --output=/path/to/stderr_stdout/folder/%x_%j.out
#SBATCH --error=/path/to/stderr_stdout/folder/%x_%j.err

################### SET VARIABLES ######################################
METADATA=metadata.R
FEATURECOUNTS=featurecounts.R

################## PRINT JOB INFO ######################################

printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

################## MODULE LOAD #########################################
#module load R/3.6.1

################## RUN SCRIPT ##########################################
Rscript $METADATA
Rscript $FEATURECOUNTS
