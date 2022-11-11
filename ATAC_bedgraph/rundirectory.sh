indir=/Shares/down/heatshock/analysis_June2022/outfiles/ATAC/ChipFlow_human/mapped/mergedbams/
outdir=/Shares/down/heatshock/analysis_June2022/outfiles/ATAC/ChipFlow_human/mapped/mergedbedgraphs/
chrom_sizes=/scratch/Shares/dowell/genomes/hg38/hg38.chrom.sizes

mkdir -p $outdir

for pathandfilename in `ls ${indir}*sorted.bam`; do
rootname=`basename $pathandfilename .sorted.bam`
echo $rootname
sbatch --export=bamfile=$pathandfilename,rootname=$rootname,chrom_sizes=$chrom_sizes,outdir=$outdir bedtoolsgeneomvecov.sh
done

