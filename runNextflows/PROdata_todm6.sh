#!/bin/bash 
#SBATCH --job-name=nextflow # Job name
#SBATCH -p long
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=1     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH --mem=8gb # Memory limit
#SBATCH --time=96:00:00 # Time limit hrs:min:sec
#SBATCH --output=/Shares/down/PRO/temp/e_and_o/nextflownascent/n.%j.out # Standard output
#SBATCH --error=/Shares/down/PRO/temp/e_and_o/nextflownascent/n.%j.err # Standard error log

nextflowdir=/Shares/down/heatshock/analysis_June2022/scripts/Nextflows/Nascent-Flow/
mkdir -p /Shares/down/PRO/temp2/
mkdir -p /Shares/down/heatshock/analysis_June2022/outfiles/PRO/Nascentflow_dros/

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


nextflow run ${nextflowdir}main.nf -profile dm6 --workdir '/Shares/down/PRO/temp2/' --outdir '/Shares/down/heatshock/analysis_June2022/outfiles/PRO/Nascentflow_dros/' --email mary.a.allen@colorado.edu --singleEnd --fastqs '/Shares/down/heatshock/analysis_June2022/fastq/PROseq_HS_lymphoblastiod/*.fastq.gz' --flip
