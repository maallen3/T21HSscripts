#!/bin/bash
#SBATCH --job-name=bidirflow                                # Job name
#SBATCH --mail-type=ALL                          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu                 # Where to send mail
#SBATCH --nodes=1                                        # Number of cores job will run on
#SBATCH --ntasks=1                                       # Number of CPU (processers, tasks)
#SBATCH --time=72:10:00                                  # Time limit hrs:min:sec
#SBATCH --partition long                                # Job queue
#SBATCH --mem=4gb                                        # Memory limit
#SBATCH --output=/scratch/Users/allenma/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/allenma/eofiles/%x_%j.err


module load samtools/1.8
module load bedtools/2.28.0
module load openmpi/1.6.4
module load gcc/7.1.0
module load python/3.6.3
module load R/3.6.1


bidirflowpath=/Shares/down/heatshock/analysis_June2022/scripts/Bidirectional-Flow/

nextflow run ${bidirflowpath}main.nf -profile hg38 \
 --crams "/Shares/down/heatshock/analysis_June2022/outfiles/PRO/Nascentflow_human/mapped/crams/*.sorted.cram" \
 --workdir "/Shares/down/PRO/temp2/" \
 --singleEnd \
 --tfit_split_model \
 --savebidirs \
 --outdir "/Shares/down/heatshock/analysis_June2022/outfiles/PRO/Nascentflow_human/mapped/bidir/"



