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
mkdir -p /Shares/down/heatshock/analysis_June2022/outfiles/ATAC/ChipFlow_human/mapped/mergedmappedfiles/
source /Users/allenma/Nascent-Flowvenv/bin/activate 


module load samtools/1.8
module load bedtools/2.28.0
module load igvtools/2.3.75
module load python/3.6.3



#nextflow run main.nf -profile hg38 --crams '/Shares/dbnascent/Allen2014global/crams/*.sorted.cram' --workdir '/scratch/Users/joca4543/191120_Cardiello_Dowell-992/results/2020_01_15/temp' --outdir '/scratch/Users/joca4543/191120_Cardiello_Dowell-992/results/2020_01_15/viz/' --saveall


nextflow run ${nextflowdir}main.nf -profile hg38 --bams '/Shares/down/heatshock/analysis_June2022/outfiles/ATAC/ChipFlow_human/mapped/mergedbams/*.sorted.bam' --workdir '/scratch/Users/allenma/temp/' --outdir '/Shares/down/heatshock/analysis_June2022/outfiles/ATAC/ChipFlow_human/mapped/mergedmappedfiles/' --saveall --resume

