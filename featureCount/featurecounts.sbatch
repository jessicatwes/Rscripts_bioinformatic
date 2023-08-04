#!/bin/bash
#SBATCH --job-name=featurecounts			 # Job name
#SBATCH --mail-type=FAIL,END				 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jewe1055@colorado.edu		 # Where to send mail
#SBATCH --nodes=1					 # Number of cores job will run on
#SBATCH --ntasks=8					 # Number of CPU (processers, tasks)
#SBATCH --time=8:00:00					 # Time limit hrs:min:sec
#SBATCH --partition short				 # Job queue
#SBATCH --mem=4gb					 # Memory limit
#SBATCH --output=/scratch/Users/jewe1055/exp_files/exp99/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/jewe1055/exp_files/exp99/eofiles/%x_%j.err

################### SET VARIABLES ######################################
INDIR='/scratch/Users/jewe1055/exp_files/exp99/scripts'
METADATA=${INDIR}/metadata.R
FEATURECOUNTS=${INDIR}/featurecounts.R

################## PRINT JOB INFO ######################################

#printf "Sample ID: $ROOTNAME"
#printf "\nDirectory: $PROJECT"
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
