
indir=/Shares/down/heatshock/analysis_June2022/outfiles/ATAC/ChipFlow_dros/mapped/bams/
outdir=/Shares/down/heatshock/analysis_June2022/outfiles/ATAC/ChipFlow_dros/mapped/mergedbams/


Vector[0]=Eric37
Vector[1]=Eric42
Vector[2]=Ethan37
Vector[3]=Ethan42


for index in $(seq 0 3)
do
rootname=${Vector[$index]}

sbatch --export=rootname=$rootname,indir=$indir,outdir=$outdir makemergeandlink.sbatch

done 
