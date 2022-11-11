library(DESeq2)
library(dplyr)

cutoff=0.1
metadatafile="/Shares/down/heatshock/analysis_June2022/scripts/T21HSscripts/metadata/PRObams_meta.txt"
coveragetablefile="/Shares/down/heatshock/analysis_June2022/outfiles/PRO/bidir_merge_diff/PRO_tfit_mumerge_peaks_165611.coverage.csv"
outdir="/Shares/down/heatshock/analysis_June2022/outfiles/PRO/bidir_merge_diff/"
indir="/Shares/down/heatshock/analysis_June2022/outfiles/PRO/bidir_merge_diff/"
annot_file = paste(indir, "PRO_seq_MUMERGE.saf", sep="")


#read the metadata
metadata <- read.csv(metadatafile, header=TRUE, sep="\t")
metadata$treatment <- relevel(metadata$tempature , "thirtyseven")
head(metadata)

safdf = read.table(annot_file, sep=",", header=TRUE)
head(safdf)
tail(safdf)

coveragetable <- read.csv(coveragetablefile, row.names=1)


#run Deseq on ATAC over peaks
ddspeaks <- DESeqDataSetFromMatrix(countData =  coveragetable, colData = metadata, design = ~genotype + treatment+genotype:treatment) #use this for all samples
dim(ddspeaks)
keep <- rowSums(counts(ddspeaks)) >= 30
ddspeaks <- ddspeaks[keep,]
dim(ddspeaks)
DEddspeaks <- DESeq(ddspeaks)
normcounts <- as.data.frame(counts(DEddspeaks, normalize=TRUE))
write.csv(as.data.frame(normcounts), file=paste(outdir,"normalizedcounts.csv",sep=""))

sample2="thirtyseven"
sample1="fourtytwo"
contrastcol = "treatment"

#get differential expresstion of regions after heat shock in D21
respeaksD21=results(DEddspeaks, contrast=c(contrastcol,sample1,sample2))
summary(respeaksD21, alpha=cutoff)
DESeq2::plotMA(respeaksD21, main="D21 differential peaks", alpha=cutoff)
write.csv(as.data.frame(respeaksD21), file=paste(outdir,"D21res.csv",sep=""))

#get differential expresstion of regions after heat shock in T21
respeaksT21=results(DEddspeaks, contrast=list( c("treatment_fourtytwo_vs_thirtyseven", "genotypeT21.treatmentfourtytwo" )))
summary(respeaksT21)
DESeq2::plotMA(respeaksT21, main="T21 differential peaks")
write.csv(as.data.frame(respeaksT21), file=paste(outdir,"T21res.csv",sep=""))

#yank only the regions different in D21 and put them in a bed files so I can look at them in igv
D21resSig_peak <- subset(respeaksD21, padj < cutoff)
D21resSig_peak$GeneID <- rownames(D21resSig_peak)
D21resSig_peak_up <- subset(D21resSig_peak, log2FoldChange > 0)
D21resSig_peak_down <- subset(D21resSig_peak, log2FoldChange < 0)


#make file bed file of upregulate regions to look at them in IGV
beddf <- merge(as.data.frame(D21resSig_peak_up),safdf,by="GeneID")
beddf <- beddf %>% dplyr::select(Chr, Start, End, GeneID, pvalue, Strand)
outfilename=paste(outdir, "D21uppeak_peakonly.unsorted.bed", sep="")
write.table(beddf, file = outfilename, sep="\t",  quote = FALSE, row.names = FALSE, col.names = FALSE)

beddf <- merge(as.data.frame(D21resSig_peak_down),safdf,by="GeneID")
beddf <- beddf %>% dplyr::select(Chr, Start, End, GeneID, pvalue, Strand)
outfilename=paste(outdir, "D21downpeak_peakonly.unsorted.bed", sep="")
write.table(beddf, file = outfilename, sep="\t",  quote = FALSE, row.names = FALSE, col.names = FALSE)

#make rank file for TFEA
respeaksD21df<- as.data.frame(respeaksD21)
respeaksD21df$GeneID <- rownames(respeaksD21df)
beddf <- merge(as.data.frame(respeaksD21df),safdf,by="GeneID")
beddf <- tibble(gene = beddf$GeneID, chrom = beddf$Chr,start = beddf$Start, stop=beddf$End,
                rnk = -log(beddf$pvalue) * sign(beddf$log2FoldChange)) %>%
  arrange(desc(rnk)) %>% drop_na()
beddf <- beddf %>% dplyr::select(chrom, start, stop) %>% dplyr::rename('#chrom'=chrom)
outfilename=paste(outdir, "D21hs_rank_file_for_TFEA.rnk", sep="")
write.table(beddf, file = outfilename, sep="\t",  quote = FALSE, row.names = FALSE)
beddf$len <- beddf$stop-beddf$start
beddf$center <- beddf$start+round(beddf$len/2, digits = 0)
beddf$start <- beddf$center-1500
beddf$stop <- beddf$center+1500
beddf <- beddf %>% dplyr::select('#chrom', start, stop) 
outfilename=paste(outdir, "D21hs_1500padded_rank_file_for_TFEA.rnk", sep="")
write.table(beddf, file = outfilename, sep="\t",  quote = FALSE, row.names = FALSE)


#yank only the regions different in D21 and put them in a bed files so I can look at them in igv
T21resSig_peak <- subset(respeaksT21, padj < cutoff)
T21resSig_peak$GeneID <- rownames(T21resSig_peak)
T21resSig_peak_up <- subset(T21resSig_peak, log2FoldChange > 0)
T21resSig_peak_down <- subset(T21resSig_peak, log2FoldChange < 0)


beddf <- merge(as.data.frame(T21resSig_peak_up),safdf,by="GeneID")
beddf <- beddf %>% dplyr::select(Chr, Start, End, GeneID, pvalue, Strand)
outfilename=paste(outdir, "T21uppeak_peakonly.unsorted.bed", sep="")
write.table(beddf, file = outfilename, sep="\t",  quote = FALSE, row.names = FALSE, col.names = FALSE)

beddf <- merge(as.data.frame(T21resSig_peak_down),safdf,by="GeneID")
beddf <- beddf %>% dplyr::select(Chr, Start, End, GeneID, pvalue, Strand)
outfilename=paste(outdir, "T21downpeak_peakonly.unsorted.bed", sep="")
write.table(beddf, file = outfilename, sep="\t",  quote = FALSE, row.names = FALSE, col.names = FALSE)

#make rank file for TFEA
respeaksT21df<- as.data.frame(respeaksT21)
respeaksT21df$GeneID <- rownames(respeaksT21df)
beddf <- merge(as.data.frame(respeaksT21df),safdf,by="GeneID")
beddf <- tibble(gene = beddf$GeneID, chrom = beddf$Chr,start = beddf$Start, stop=beddf$End,
                rnk = -log(beddf$pvalue) * sign(beddf$log2FoldChange)) %>%
  arrange(desc(rnk)) %>% drop_na()
beddf <- beddf %>% dplyr::select(chrom, start, stop) %>% dplyr::rename('#chrom'=chrom)
outfilename=paste(outdir, "T21hs_rank_file_for_TFEA.rnk", sep="")
write.table(beddf, file = outfilename, sep="\t",  quote = FALSE, row.names = FALSE)
beddf$len <- beddf$stop-beddf$start
beddf$center <- beddf$start+round(beddf$len/2, digits = 0)
beddf$start <- beddf$center-1500
beddf$stop <- beddf$center+1500
beddf <- beddf %>% dplyr::select('#chrom', start, stop) 
outfilename=paste(outdir, "T21hs_1500padded_rank_file_for_TFEA.rnk", sep="")
write.table(beddf, file = outfilename, sep="\t",  quote = FALSE, row.names = FALSE)
