#!/bin/bash
#SBATCH --job-name=idxstats
#SBATCH --mail-type=ALL
#SBATCH --mail-user=%u@colorado.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100GB
#SBATCH --time=23:59:59
#SBATCH --partition=short
#SBATCH --output=/scratch/Users/allenma/eofiles/%x.%j.out
#SBATCH --error=/scratch/Users/allenma/eofiles/%x.%j.err

echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host `hostname`
echo Job started at `date +"%T %a %d %b %Y"`
echo Directory is `pwd`
echo Using $SLURM_NTASKS processors, across $SLURM_NNODES nodes, with $SLURM_JOB_CPUS_PER_NODE cpus per node

module purge


bamdir=/Shares/down/heatshock/analysis_June2022/outfiles/PRO/BedgraphsandBigwigs_human/mapped/bams/
outfile=/Shares/down/heatshock/analysis_June2022/outfiles/PRO/stats/idxstats.txt

source ~/idxstatsvenv/bin/activate 

python3 readsperch.py $bamdir $outfile
