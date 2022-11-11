#!/bin/bash 
#SBATCH --job-name=reindex # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=1     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH --mem=10gb # Memory limit
#SBATCH --time=23:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Users/allenma/eofiles/reindex.%j.out # Standard output
#SBATCH --error=/scratch/Users/allenma/eofiles/reindex.%j.err # Standard error log


module load samtools/1.8


samtools index ${pathandfilename} 
