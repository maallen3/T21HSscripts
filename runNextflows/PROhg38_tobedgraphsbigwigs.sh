#!/bin/bash 
#SBATCH --job-name=nextflow # Job name
#SBATCH -p long
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=1     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH --mem=8gb # Memory limit
#SBATCH --time=96:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Users/allenma/e_and_o/nextflow.%j.out # Standard output
#SBATCH --error=/scratch/Users/allenma/e_and_o/nextflow.%j.err # Standard error log

nextflowdir=/Shares/down/heatshock/analysis_June2022/scripts/Nextflows/Downfile_pipeline/
mkdir -p /Shares/down/heatshock/analysis_June2022/outfiles/PRO/BedgraphsandBigwigs_human/
source /Shares/down/heatshock/analysis_June2022/scripts/venv_Nextflow/bin/activate

module load samtools/1.8
module load bedtools/2.28.0
module load igvtools/2.3.75
module load python/3.6.3


nextflow run ${nextflowdir}main.nf -profile hg38 --crams '/Shares/down/heatshock/analysis_June2022/outfiles/PRO/Nascentflow_human/mapped/crams/*.sorted.cram' --workdir '/scratch/Users/allenma/temp/' --outdir '/Shares/down/heatshock/analysis_June2022/outfiles/PRO/BedgraphsandBigwigs_human/' --saveall --singleEnd

