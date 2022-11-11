#!/bin/bash
#SBATCH --job-name=pyGenomeTracks				 # Job name
#SBATCH --mail-type=ALL				 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu		 # Where to send mail
#SBATCH --nodes=1					 # Number of cores job will run on
#SBATCH --ntasks=1					 # Number of CPU (processers, tasks)
#SBATCH --time=00:10:00					 # Time limit hrs:min:sec
#SBATCH --partition short				 # Job queue
#SBATCH --mem=4gb					 # Memory limit
#SBATCH --output=/scratch/Users/allenma/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/allenma/eofiles/%x_%j.err

################### SET VARIABLES ######################################
OUT='/Shares/down/heatshock/analysis_June2022/outfiles/PRO/igv_combined'
INI='/Shares/down/heatshock/analysis_June2022/scripts/T21HSscripts/gene_figures/ini'
ROOTNAME='pro_t21_hs_serpinh1'

########################################################################
################### LOAD NECESSARY MODULES #############################
module load gcc/7.2.0 
module load bedtools ##this is needed to add the genome track
source /Shares/down/heatshock/analysis_June2022/scripts/T21HSscripts/gene_figures/venv_pygenometracks/bin/activate


########################################################################
################## PRINT JOB INFO ######################################

printf "Output Dir: $OUT"
printf "\nINI File: $INI"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

#######################################################################

pyGenomeTracks \
 --tracks ${INI}/${ROOTNAME}.ini \
 --region chr11:75,561,294-75,582,876 \
 --outFileName $OUT/${ROOTNAME}.svg

printf "complete"
