#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
library("DESeq2")
library("pheatmap")
library("dplyr")

#set your out directory and the gtf for your genome
indir <- "/Shares/down/heatshock/analysis_June2022/outfiles/PRO/gene_counts_and_diffexp_removechr/body_only/"
outdir="/Shares/down/heatshock/analysis_June2022/outfiles/PRO/gene_counts_and_diffexp_removechr/body_only/"
annot_file <-"/scratch/Shares/dowell/genomes/hg38/PROparts_genenames/Nov2022_removechrnotinfasta/body.saf"
metadir <- "/Shares/down/heatshock/analysis_June2022/scripts/T21HSscripts/metadata/"
cutoff = 0.01

setwd(outdir)

#read the counts csv
fileroot = "featureCounts_body_112233"
coveragetable <- read.csv(paste(indir, fileroot, ".coverage.csv", sep=""), row.names=1)

head(coveragetable)

#read the metadata
metadata <- read.csv(paste(metadir, "PRObams_meta.txt", sep=""), header=TRUE, sep="\t")
head(metadata)
metadata$genotype <- relevel(metadata$genotype , "D21")
head(metadata)
metadata$tempature <- relevel(metadata$tempature , "thirtyseven")
head(metadata)

#make sure they are in the same order and have the same samples!!!!!
countdat <- coveragetable %>% dplyr::select(as.vector(metadata$label))
dim(countdat)

#set up the deseq object
dds <- DESeqDataSetFromMatrix(countData = countdat, colData = metadata, design = ~group)
ddshs <- DESeqDataSetFromMatrix(countData = countdat, colData = metadata, design = ~tempature)


#run deseq
DEdds <- DESeq(dds)
DEddshs <- DESeq(ddshs)

reshs <- results(DEddshs)

#Check the size factors
sizeFactors(DEdds) #this should be pretty simlilar to ratios in the millions mapped. CHECK IT!


#these two graphs show you before and after normalization of samples. In this data set that doesn't matter. In some it does.
allcounts <- as.data.frame(counts(DEdds))
allcountslong <- allcounts %>% gather(key = "sample", value = "signal")
ggplot(allcountslong, aes(x = sample, y = signal))+ geom_violin(trim = FALSE) + theme(axis.text.x = element_text(angle = 90))+scale_y_continuous(trans='log2')

normcounts <- as.data.frame(counts(DEdds, normalize=TRUE))
normcountslong <- normcounts %>% gather(key = "sample", value = "signal")
ggplot(normcountslong, aes(x = sample, y = signal))+ geom_violin(trim = FALSE) + theme(axis.text.x = element_text(angle = 90))+scale_y_continuous(trans='log2')

#plot the dispersion of the data
plotDispEsts(DEdds)

#get the results dataframes from comparing two samples 
#my samples are D21_fourtytwo D21_thirtyseven T21_fourtytwo T21_thirtyseven
#Then save the results file


outputfiles <-function(res, sample1, sample2){
  res_Sig <- subset(res, padj < cutoff)
  write.csv(res, file = paste0(outdir,sample1, "_",sample2, "_all_results.csv"))
  write.csv(res_Sig, file = paste0(outdir,sample1, "_",sample2, "_sig_results.csv"))
  #for go and enricher and gsea
  res_expressed <- res[!is.na(res$padj),]
  write.csv(rownames(res_expressed), file = paste0(outdir,sample1, "_",sample2, "T21_backgroundgenes.csv"),row.names = FALSE, quote = FALSE)
  res_expressed_sig <- subset(res_expressed, padj < cutoff)
  write.csv(rownames(res_expressed_sig), file = paste0(outdir,sample1, "_",sample2,"T21_siggenes.csv"),row.names = FALSE, quote = FALSE)
  rnkdf <- tibble(gene = rownames(res),
                  rnk = -log(res$pvalue) * sign(res$log2FoldChange)) %>%arrange(desc(rnk)) %>% drop_na()
  write.table(rnkdf, file = paste0(outdir,sample1, "_",sample2, "_deseq_res_for_gsea.rnk"),
              append = FALSE, col.names = FALSE, row.names = FALSE,
              quote = FALSE, sep = "\t")}

sample1="T21_thirtyseven"
sample2="T21_fourtytwo"
res_T21=results(DEdds, contrast=c("group",sample1,sample2))
write.csv(res_T21, paste(outdir,fileroot,"_",sample1, "_",sample2,".results.csv", sep=""))
outputfiles(res_T21, sample1, sample2)

sample1="D21_thirtyseven"
sample2="D21_fourtytwo"
res_D21=results(DEdds, contrast=c("group",sample1,sample2))
write.csv(res_D21, paste(outdir,fileroot,"_",sample1, "_",sample2,".results.csv", sep=""))
outputfiles(res_D21, sample1, sample2)

#look at the results
res<-res_T21
summary(res)
DESeq2::plotMA(res)
resSig <- subset(res, padj < cutoff)
head(res[ order( res$padj ), ])
dim(resSig)

#look at the results
res<-res_D21
summary(res)
DESeq2::plotMA(res)
resSig <- subset(res, padj < cutoff)
head(res[ order( res$padj ), ])
dim(resSig)


#plot one gene
#genename="NM_020161"
#genename=" NM_001017963_"
#plotCounts(DEdds, gene=genename, intgroup=c("group"))

plotlogfc<- function(res1, res2, name1, name2){
  res1 = as.data.frame(res1)
  res2 = as.data.frame(res2)
  df <- merge(res1, res2, by.x=0, by.y=0, suffixes=c(paste(".",name1, sep=""), paste(".",name2, sep="")))
  }

df <- plotlogfc(res_D21, res_T21, "D21", "T21")
df_sig <- df %>% filter(padj.T21<cutoff | padj.D21<cutoff)
df_sig <- df_sig %>%  mutate(max_padj = pmax(padj.T21,padj.D21))
df_sig <- df_sig %>%  mutate(min_padj = pmin(padj.T21,padj.D21))
df_sig <- df_sig %>% mutate(cat = case_when(max_padj<cutoff~"both", padj.T21 <cutoff ~"T21", padj.D21 <cutoff ~"D21"))

df <- df %>%  mutate(max_padj = pmax(padj.T21,padj.D21))
df <- df %>%  mutate(min_padj = pmin(padj.T21,padj.D21))
df <- df %>% mutate(cat = case_when(max_padj <cutoff~"both", padj.T21 <cutoff ~"T21", padj.D21 <cutoff ~"D21", max_padj >cutoff~"N.S"))


df_sig %>% group_by(cat) %>% summarize(count=n()) 
df_sig_D21_only_up <- df_sig %>% filter(cat=="D21") %>% filter(log2FoldChange.D21>0)
df_sig_D21_only_down <- df_sig %>% filter(cat=="D21")%>% filter(log2FoldChange.D21<0)
df_sig_T21_only_up <- df_sig %>% filter(cat=="T21") %>% filter(log2FoldChange.T21>0)
df_sig_T21_only_down <- df_sig %>% filter(cat=="T21") %>% filter(log2FoldChange.T21<0)

head(df_sig_D21_only_up[ order( df_sig_D21_only_up$padj.D21 ), ])

df %>% group_by(cat) %>% summarize(count=n()) 


ggplot(df_sig, aes(x=padj.T21, y=padj.D21, color=cat))+geom_point()
ggplot(df_sig, aes(x=log10(padj.T21), y=log10(padj.D21), color=cat))+geom_point()

ggplot(df_sig, aes(x=log10(padj.T21), y=log10(padj.D21), color=cat))+geom_point()+xlim(-10,0)+ylim(-10,0)

ggplot(df_sig, aes(x=log2FoldChange.D21, y=log2FoldChange.T21, color=cat))+geom_point()+geom_abline(slope=1, intercept=0)

ggplot(df_sig, aes(x=log2FoldChange.D21, y=log2FoldChange.T21, color=cat))+geom_point()+geom_abline(slope=1, intercept=0)+facet_wrap(~cat)


hist(df_sig$padj.T21, breaks=100)
hist(df_sig$padj.D21, breaks=100)

ggplot(df, aes(x=baseMean.D21, y=log2FoldChange.D21, color=cat))+geom_point()+scale_x_continuous(trans="log10")+facet_wrap(~cat)
ggplot(df, aes(x=baseMean.T21, y=log2FoldChange.T21, color=cat))+geom_point()+scale_x_continuous(trans="log10")+facet_wrap(~cat)

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}



normcounts_group <- normcounts[rownames(normcounts) %in% noquote(row.names(resSig)),]
normcounts_group <- as.matrix(normcounts_group)
normcounts_group_rowscaled <- t(apply(normcounts_group, 1, cal_z_score))
pheatmap(normcounts_group_rowscaled, cluster_cols = F)

normcounts_group <- normcounts[rownames(normcounts) %in% noquote(df_sig$Row.names),]
normcounts_group <- as.matrix(normcounts_group)
normcounts_group_rowscaled <- t(apply(normcounts_group, 1, cal_z_score))
pheatmap(normcounts_group_rowscaled, cluster_cols = F)

normcounts_group <- normcounts[rownames(normcounts) %in% noquote(df_sig_D21_only_up$Row.names),]
normcounts_group <- as.matrix(normcounts_group)
normcounts_group_rowscaled <- t(apply(normcounts_group, 1, cal_z_score))
pheatmap(normcounts_group_rowscaled, cluster_cols = F)

normcounts_group <- normcounts[rownames(normcounts) %in% noquote(df_sig_D21_only_down$Row.names),]
normcounts_group <- as.matrix(normcounts_group)
normcounts_group_rowscaled <- t(apply(normcounts_group, 1, cal_z_score))
pheatmap(normcounts_group_rowscaled, cluster_cols = F)

normcounts_group <- normcounts[rownames(normcounts) %in% noquote(df_sig_T21_only_up$Row.names),]
normcounts_group <- as.matrix(normcounts_group)
normcounts_group_rowscaled <- t(apply(normcounts_group, 1, cal_z_score))
pheatmap(normcounts_group_rowscaled, cluster_cols = F)

normcounts_group <- normcounts[rownames(normcounts) %in% noquote(df_sig_T21_only_down$Row.names),]
normcounts_group <- as.matrix(normcounts_group)
normcounts_group_rowscaled <- t(apply(normcounts_group, 1, cal_z_score))
pheatmap(normcounts_group_rowscaled, cluster_cols = F)




onegene = "NM_005526_HSF1"
plotCounts(dds, onegene, intgroup=c("genotype", "tempature"))
plotCounts(dds, onegene, intgroup=c("replicate", "genotype", "tempature"))
