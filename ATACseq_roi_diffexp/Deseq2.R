
library(DESeq2)
library(dplyr)

cutoff=0.1

metadatafile="/Shares/down/heatshock/analysis_June2022/scripts/T21HSscripts/metadata/ATACbams_meta.txt"
indir="/Shares/down/heatshock/analysis_June2022/outfiles/ATAC/peaks_counts_and_diff/"
coveragetablefile=paste(indir, "featureCounts_open_122800.coverage.csv", sep="")
mumergedir = "/Shares/down/heatshock/analysis_June2022/outfiles/ATAC/peak_merge/"
annot_file =paste(mumergedir, "ATACHMMRATACmerge_MUMERGE.saf", sep="")
outdir="/Shares/down/heatshock/analysis_June2022/outfiles/ATAC/peaks_counts_and_diff/"

#read the metadata
metadata <- read.csv(metadatafile, header=TRUE, sep="\t")
head(metadata)

metadata$genotype <- relevel(metadata$genotype , "D21")
metadata$tempature <- relevel(metadata$tempature , "thirtyseven")
head(metadata)


safdf = read.table(annot_file, sep="\t", header=TRUE)
head(safdf)
tail(safdf)


coveragetable  <- read.csv(coveragetablefile, row.names=1)
# this rearagnes the rows so they are in the same order as the metadata
countdat <- coveragetable %>% dplyr::select(as.vector(metadata$label))



dds <- DESeqDataSetFromMatrix(countData =  countdat, colData = metadata, design = ~group) #use this for all samples

keep <- rowSums(counts(dds)) >= 30
dds <- dds[keep,]
DEdds <- DESeq(dds)
normcounts <- as.data.frame(counts(DEdds, normalize=TRUE))

sf <- sizeFactors(DEdds)
write.csv(as.data.frame(sf), file=paste(outdir,"sizefactors.csv",sep=""),quote=FALSE)


sample2="D21_thirtyseven"
sample1="D21_fourtytwo"
contrastcol = "group"

respeaksD21=results(DEdds, contrast=c(contrastcol,sample1,sample2))
summary(respeaksD21, alpha=cutoff)
DESeq2::plotMA(respeaksD21, main="D21 differential peaks", alpha=cutoff)
write.csv(as.data.frame(respeaksD21), file=paste(outdir,"D21res.csv",sep=""))

sample2="T21_thirtyseven"
sample1="T21_fourtytwo"
contrastcol = "group"

respeaksT21=results(DEdds, contrast=c(contrastcol,sample1,sample2))
summary(respeaksT21, alpha=cutoff)
DESeq2::plotMA(respeaksT21, main="T21 differential peaks", alpha=cutoff)
write.csv(as.data.frame(respeaksT21), file=paste(outdir,"T21res.csv",sep=""))

D21resSig_peak <- subset(respeaksD21, padj < cutoff)
D21resSig_peak$GeneID <- rownames(D21resSig_peak)
D21resSig_peak_up <- subset(D21resSig_peak, log2FoldChange > 0)
D21resSig_peak_down <- subset(D21resSig_peak, log2FoldChange < 0)

head(D21resSig_peak_up[ order( D21resSig_peak_up$padj ), ])
head(D21resSig_peak_down[ order( D21resSig_peak_down$padj ), ])

beddf <- merge(as.data.frame(D21resSig_peak_up),safdf,by="GeneID")
beddf <- beddf %>% dplyr::select(Chr, Start, End, GeneID, pvalue, Strand)
outfilename=paste(outdir, "D21uppeak_peakonly.unsorted.bed", sep="")
write.table(beddf, file = outfilename, sep="\t",  quote = FALSE, row.names = FALSE, col.names = FALSE)

beddf <- merge(as.data.frame(D21resSig_peak_down),safdf,by="GeneID")
beddf <- beddf %>% dplyr::select(Chr, Start, End, GeneID, pvalue, Strand)
outfilename=paste(outdir, "D21downpeak_peakonly.unsorted.bed", sep="")
write.table(beddf, file = outfilename, sep="\t",  quote = FALSE, row.names = FALSE, col.names = FALSE)



T21resSig_peak <- subset(respeaksT21, padj < cutoff)
T21resSig_peak$GeneID <- rownames(T21resSig_peak)
T21resSig_peak_up <- subset(T21resSig_peak, log2FoldChange > 0)
T21resSig_peak_down <- subset(T21resSig_peak, log2FoldChange < 0)

head(T21resSig_peak_up[ order( T21resSig_peak_up$padj ), ])
head(T21resSig_peak_down[ order( T21resSig_peak_down$padj ), ])

beddf <- merge(as.data.frame(T21resSig_peak_up),safdf,by="GeneID")
beddf <- beddf %>% dplyr::select(Chr, Start, End, GeneID, pvalue, Strand)
outfilename=paste(outdir, "T21uppeak_peakonly.unsorted.bed", sep="")
write.table(beddf, file = outfilename, sep="\t",  quote = FALSE, row.names = FALSE, col.names = FALSE)

beddf <- merge(as.data.frame(T21resSig_peak_down),safdf,by="GeneID")
beddf <- beddf %>% dplyr::select(Chr, Start, End, GeneID, pvalue, Strand)
outfilename=paste(outdir, "T21downpeak_peakonly.unsorted.bed", sep="")
write.table(beddf, file = outfilename, sep="\t",  quote = FALSE, row.names = FALSE, col.names = FALSE)




