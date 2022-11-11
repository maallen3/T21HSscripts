#!/bin/bash
#SBATCH --job-name=mumerge                     # Job name
#SBATCH --mail-type=ALL                           # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu                       # Where to send mail
#SBATCH --nodes=1                               # Number of nodes requested
#SBATCH --ntasks=1                              # Number of CPUs (processor cores/tasks)
#SBATCH --mem=8gb                               # Memory limit
#SBATCH --time=23:30:00                         # Time limit hrs:min:sec
#SBATCH --partition=short                       # Partition/queue requested on server
#SBATCH --output=/scratch/Users/allenma/e_and_o/%x.%j.out
#SBATCH --error=/scratch/Users/allenma/e_and_o/%x.%j.err

module load bedtools/2.25.0

source /Shares/down/heatshock/analysis_June2022/scripts/vir_TFEA/bin/activate
mumergedir=/Shares/down/heatshock/analysis_June2022/scripts/mumerge/mumerge/
sample_information=/Shares/down/heatshock/analysis_June2022/scripts/T21HSscripts/metadata/ATACHMMRbeds.txt
out_path_project_id=/Shares/down/heatshock/analysis_June2022/outfiles/ATAC/peak_merge/ATACHMMRATACmerge

python ${mumergedir}mumerge.py -i $sample_information -o $out_path_project_id
