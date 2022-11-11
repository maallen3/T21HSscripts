

beddir=/Shares/down/heatshock/analysis_June2022/outfiles/ATAC/peaks_counts_and_diff/
genefile=/scratch/Shares/dowell/genomes/hg38/hg38_refseq_genenames.bed
#need unsortedbedfile, peakfile, genesnearpeaksfile
outdir=$beddir

for pathandfilename in `ls ${beddir}*.unsorted.bed`; do
rootname=`basename $pathandfilename .unsorted.bed`
peakfile=${outdir}${rootname}.sorted.bed
mkdir -p ${outdir}/nearbygenes/
genesnearpeaksfile=${outdir}/nearbygenes/genesnear_${rootname}.bed
echo $rootname
sbatch --export=unsortedbedfile=$pathandfilename,peakfile=$peakfile,genesnearpeaksfile=$genesnearpeaksfile,genefile=$genefile /Shares/down/heatshock/analysis_June2022/scripts/T21HSscripts/ATACseq_roi_diffexp/diffpeaksneargenes.sh
done

