#!/bin/bash
#SBATCH --job-name=tfea
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

export PATH=~:$PATH
export PATH=~/.local/bin:$PATH

module load python/3.7.4
module load samtools/1.3.1
module load bedtools/2.25.0
module load meme/5.0.3
module load samtools/1.3.1
module load gcc/7.1.0
module load R/4.2.0



source /Shares/down/heatshock/analysis_June2022/scripts/vir_TFEA/bin/activate
rank_dir=/Shares/down/heatshock/analysis_June2022/outfiles/PRO/bidir_merge_diff/
rank_file=${rank_dir}T21hs_rank_file_for_TFEA.rnk
outdir=/Shares/down/heatshock/analysis_June2022/outfiles/PRO/TFEA_output/HOCOMOCO10/Ethan_hs_via_rank/
TFEAdir=/Shares/down/heatshock/analysis_June2022/scripts/TFEA
fastadir=/scratch/Shares/dowell/genomes/hg38/
HOCOMOCOdir=/scratch/Shares/dowell/motifs/HOCOMOCODatabaseFIMO/

mkdir -p $outdir

python3 /Shares/down/heatshock/analysis_June2022/scripts/vir_TFEA/lib/python3.7/site-packages/tfea-1.1.4-py3.7.egg/TFEA \
--output ${outdir} \
--ranked_file $rank_file \
--label1 T21_37 \
--label2 T21_42 \
--genomefasta ${fastadir}hg38.fa \
--fimo_motifs ${HOCOMOCOdir}HOCOMOCOv10_HUMAN_mono_meme_format.meme \
--cpus 16 \
--mem 100gb \
--output_type html \



