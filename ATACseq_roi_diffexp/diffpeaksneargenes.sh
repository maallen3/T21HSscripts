#!/bin/bash 
#SBATCH --job-name=peaksneargenes # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=1     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH --mem=10gb # Memory limit
#SBATCH --time=23:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Users/allenma/e_and_o/peaksneargenes.%j.out # Standard output
#SBATCH --error=/scratch/Users/allenma/e_and_o/peaksneargenes.%j.err # Standard error log


module load gcc/7.2.0 
module load bedtools/2.28.0 
source /Shares/down/heatshock/analysis_June2022/scripts/vir_TFEA/bin/activate

bedtools sort -i $unsortedbedfile >$peakfile 
bedtools window -w 25000 -u -a $genefile -b $peakfile >$genesnearpeaksfile
python yankgenename.py $genesnearpeaksfile
python yankgenename.py $genefile

