#/bin/bash
#SBATCH --job-name=mergeandsort                     # Job name
#SBATCH --mail-type=ALL                           # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu                       # Where to send mail
#SBATCH --nodes=1                               # Number of nodes requested
#SBATCH --ntasks=32                              # Number of CPUs (processor cores/tasks)
#SBATCH --mem=256gb                               # Memory limit
#SBATCH --time=23:30:00                         # Time limit hrs:min:sec
#SBATCH --partition=short                       # Partition/queue requested on server
#SBATCH --output=/scratch/Users/allenma/e_and_o/%x.%j.out
#SBATCH --error=/scratch/Users/allenma/e_and_o/%x.%j.err 

### Displays the job context
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host `hostname`
echo Job started at `date +"%T %a %d %b %Y"`
echo Directory is `pwd`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

pwd; hostname; date

module load samtools/1.8 

ln -sf ${indir}${rootname}_rep1.sorted.bam ${outdir}${rootname}_rep1.sorted.bam
ln -sf ${indir}${rootname}_rep1.sorted.bam.bai ${outdir}${rootname}_rep1.sorted.bam.bai
#mergeing the technical replicates
samtools merge -@ 32 ${outdir}${rootname}_rep2.bam ${indir}${rootname}_rep2a.sorted.bam ${indir}${rootname}_rep2b.sorted.bam
samtools sort ${outdir}${rootname}_rep2.bam -m 8G -@ 32 -O bam -o ${outdir}${rootname}_rep2.sorted.bam
samtools index  ${outdir}${rootname}_rep2.sorted.bam
rm ${outdir}${rootname}_rep2.bam



