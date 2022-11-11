#!/bin/bash 
#SBATCH --job-name=nextflow # Job name
#SBATCH -p long
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=1     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH --mem=8gb # Memory limit
#SBATCH --time=96:00:00 # Time limit hrs:min:sec
#SBATCH --output=/Shares/down/ATAC/temp/e_and_o/nextflowCHIP.%j.out # Standard output
#SBATCH --error=/Shares/down/ATAC/temp/e_and_o/nextflowCHIP.%j.err # Standard error log

nexflowdir=/Shares/down/heatshock/analysis_June2022/scripts/Nextflows/ChIP-Flow/
mkdir -p /Shares/down/ATAC/temp2/
mkdir -p /Shares/down/heatshock/analysis_June2022/outfiles/ATAC/ChipFlow_dros/

source /Shares/down/heatshock/analysis_June2022/scripts/venv_Nextflow/bin/activate

module load sra/2.9.2
module load bbmap/38.05
module load fastqc/0.11.8
module load hisat2/2.1.0
module load samtools/1.8
module load preseq/2.0.3
module load python/3.6.3/rseqc
module load bedtools/2.28.0
module load igvtools/2.3.75
module load mpich/3.2.1
module load openmpi/1.6.4
module load gcc/7.2.0


nextflow run ${nexflowdir}main.nf -profile slurm_dro --workdir '/Shares/down/ATAC/temp2/' --outdir '/Shares/down/heatshock/analysis_June2022/outfiles/ATAC/ChipFlow_dros/' --email mary.a.allen@colorado.edu --fastqs '/Shares/down/heatshock/analysis_June2022/fastq/ATAC_HS_lymphoblastiod/*{R1,R2}*.fastq.gz'
