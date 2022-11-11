indir=/Shares/down/heatshock/analysis_June2022/outfiles/ATAC/ChipFlow_dros/mapped/mergedbams/
outdir=/Shares/down/heatshock/analysis_June2022/outfiles/ATAC/peaks/HMMRATACpeaks_dm6/
HMMRATACdir=/Shares/down/heatshock/analysis_June2022/scripts/HMMRATAC/

mkdir -p $outdir


for pathandfilename in `ls ${indir}*sorted.bam`; do
rootname=`basename $pathandfilename .sorted.bam`
echo $indir
echo $rootname
sbatch --export=indir=$indir,rootname=$rootname,outdir=$outdir,HMMRATACdir=$HMMRATACdir HMMRATACandsort_python.sh
done

echo results wiil be in $outdir
