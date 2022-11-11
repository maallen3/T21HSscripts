#indir=/Shares/down/heatshock/analysis_June2022/outfiles/ATAC/ChipFlow_dros/mapped/bams/
#indir=/Shares/down/heatshock/analysis_June2022/outfiles/PRO/BedgraphsandBigwigs_human/mapped/bams/
indir=/Shares/down/heatshock/analysis_June2022/outfiles/PRO/BedgraphsandBigwigs_dros/mapped/bams/

for pathandfilename in `ls ${indir}*sorted.bam`; do
sbatch --export=pathandfilename=$pathandfilename samtoolsindex.sh
done

