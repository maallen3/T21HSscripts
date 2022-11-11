from __future__ import division
import pysam
import pandas
import glob
import sys
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

allchrs = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', "chrM"]

allchrs1000g = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', "MT"]


def idxcountsperchr(bamfile):
	lensdic={}
	obamfile = pysam.AlignmentFile(bamfile, "rb")
	h= obamfile.header
	for stuff in  h['SQ']:
		lenofchr = stuff['LN']
		chrname = stuff['SN']
		lensdic[chrname] = lenofchr
	stats =obamfile.get_index_statistics() 
	totalmappedreads = sum([x.total for x in stats])
	mm = totalmappedreads/1000000
	chrnames = [[x.contig,x.mapped,x.total,lensdic[x.contig],(x.mapped/(lensdic[x.contig]*mm))] for x in stats]
	return chrnames	




def idxcountsdir_total(dirname, writefile):
	dfs = []
	for bamfile in sorted(glob.glob(dirname+"*.bam")):
		print(bamfile)
		chrnames= idxcountsperchr(bamfile)
		rootname = bamfile.split("/")[-1]
		#rootname = rootname.split("_")[0]
		labels = ["chrname", rootname+"_mapped", rootname+"_total", rootname+"_lenchr", rootname+"_fpkm"]
		df = pandas.DataFrame.from_records(chrnames,columns=labels)
		df.set_index('chrname', inplace=True)	
		dfs.append(df)
	print (dfs)
	result = pandas.concat(dfs, axis=1)
	result.to_csv(writefile+".fullinfo.txt")
	final = result.filter(regex='total', axis=1)
	final.to_csv(writefile)
	print (final)
	uallchrs = [str(c) for c in allchrs if c!="chrM"]+[str(c) for c in allchrs1000g if c!="MT"]
	finallimitchr= final[final.index.map(lambda x: x in uallchrs)]
	fig, axs = plt.subplots(ncols=2)
	sns.heatmap(finallimitchr, ax = axs[0])
	sns.heatmap(finallimitchr, ax = axs[1],norm=LogNorm(vmin=finallimitchr.min(), vmax=finallimitchr.max()))
	axs[0].set_title('Linear norm colorbar, seaborn')
	axs[1].set_title('Log norm colorbar, seaborn')	
	plt.subplots_adjust(bottom=0.4)
	#plt.tight_layout()
	plt.savefig(writefile+".png")	
	





def graphfpkmfc(idxstatfpkmfile="/scratch/Shares/dowell/Down/genome/paperfigs042518/graphs/chrfpkmourpeople.txt", foldchangedenom="WandaGenome.bam_fpkm",writefile="/scratch/Shares/dowell/Down/genome/paperfigs042518/graphs/chrfpkmourpeoplefc"):
	df = pandas.read_csv(idxstatfpkmfile,index_col=0)
	df = df.iloc[:,1:].div(df[foldchangedenom], axis=0)
	final = df
	uallchrs = [str(c) for c in allchrs if c!="chrM"]+[str(c) for c in allchrs1000g if c!="MT"]
	finallimitchr= final[final.index.map(lambda x: x in uallchrs)]
	fig, axs = plt.subplots(ncols=2)
	sns.heatmap(finallimitchr, ax = axs[0])
	sns.heatmap(finallimitchr, ax = axs[1],norm=LogNorm(vmin=finallimitchr.min(), vmax=finallimitchr.max()))
	axs[0].set_title('Linear norm colorbar, seaborn')
	axs[1].set_title('Log norm colorbar, seaborn')
	plt.tight_layout()
	plt.savefig(writefile+".png")


if __name__=="__main__":
	indir = sys.argv[1]
	writefile = sys.argv[2]
	idxcountsdir_total(indir, writefile)

