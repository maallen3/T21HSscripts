
library(DESeq2)
library(dplyr)
library(DEGreport)

cutoff=0.1
metadatafile="/Shares/down/heatshock/analysis_June2022/scripts/T21HSscripts/metadata/ATACbams_meta.txt"
indir="/Shares/down/heatshock/analysis_June2022/outfiles/ATAC/peaks_counts_and_diff/"
coveragetablefile=paste(indir, "featureCounts_open_122800.coverage.csv", sep="")
annot_file =paste(mumergedir, "ATACHMMRATACmerge_MUMERGE.saf", sep="")
outdir="/Shares/down/heatshock/analysis_June2022/outfiles/ATAC/peaks_counts_and_diff_LRT/"

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


dds <- DESeqDataSetFromMatrix(countData =  countdat, colData = metadata, design = ~genotype+tempature+genotype:tempature) #use this for all samples



dds1 <- DESeqDataSetFromMatrix(countData =  countdat, colData = metadata, design = ~genotype+tempature+genotype:tempature) #use this for all samples
dds2 <- DESeqDataSetFromMatrix(countData =  countdat, colData = metadata, design = ~genotype+tempature+genotype:tempature) #use this for all samples
dds3 <- DESeqDataSetFromMatrix(countData =  countdat, colData = metadata, design = ~genotype+tempature+genotype:tempature) #use this for all samples
dds4 <- DESeqDataSetFromMatrix(countData =  countdat, colData = metadata, design = ~genotype+tempature+genotype:tempature) #use this for all samples


dds_lrt_1 <- DESeq(dds1, test="LRT", reduced = ~ 1)
dds_lrt_t <- DESeq(dds2, test="LRT", reduced = ~tempature)
dds_lrt_g <- DESeq(dds3, test="LRT", reduced = ~genotype)
dds_lrt_gt <- DESeq(dds4, test="LRT", reduced = ~genotype+tempature)


res_1 <- results(dds_lrt_1)
resOrdered_1 <- res_1[order(res_1$pvalue),]

res_t <- results(dds_lrt_t)
resOrdered_t <- res_t[order(res_t$pvalue),]

res_g <- results(dds_lrt_g)
resOrdered_g <- res_g[order(res_g$pvalue),]

res_gt <- results(dds_lrt_gt)
resOrdered_gt <- res_gt[order(res_gt$pvalue),]
hist(resOrdered_gt$pvalue, breaks=100)


res_LRT <- res_1


padj.cutoff <-0.01
sig_res_LRT <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff)

# Get sig gene lists
sigLRT_genes <- sig_res_LRT %>% 
  pull(gene)

length(sigLRT_genes)

#clustering takes a while so we are only going to cluster on the top 1000 genes in class
#make sure to use all genes when running this for real
clustering_sig_genes <- sig_res_LRT %>%
  arrange(padj) %>%
  head(n=1000)

# Obtain rlog values for those significant genes
rld <- rlog(dds_lrt_1, blind=TRUE)
rld_mat <- assay(rld)
cluster_rlog <- rld_mat[sig_res_LRT$gene, ]

#make a metadata frame for this anaysis... different becuase the index must be the column names for the counts
meta <-as.data.frame(cbind(paste(dds_lrt_1$label), paste(dds_lrt_1$genotype),paste(dds_lrt_1$tempature),paste(DEdds$group), paste(DEdds$replicate)))
colnames(meta)<-c("label","genotype",  "tempature", "group", "replicate")
row.names(meta)<-meta$label
meta$tempature <-relevel(meta$tempature , "thirtyseven")

#dev.off()
#jpeg(paste0(outdir,outfilename,"/",'clusters.jpg', sep=""))
#clusters <- degPatterns(cluster_rlog, metadata = meta, time = "group", col = NULL)
clusters <- degPatterns(cluster_rlog, metadata = meta, time = "tempature", col = "genotype",reduce=TRUE)

#dev.off()

#get more info on clusters
cluster_groups<- clusters$df
for (i in unique(cluster_groups$cluster)){
  group <- clusters$df %>% filter(cluster == i)
  normcounts_group <- normcounts[rownames(normcounts) %in% noquote(group$genes),]
  normcounts_group <- as.matrix(normcounts_group)
  title = paste0("cluster", i, sep=" ")
  heatmap(normcounts_group,Colv = NA, main=title)}


maindf <- merge(as.data.frame(res_LRT), cluster_groups, by=0)

cluster1 = maindf %>% filter(cluster==1) %>% arrange(padj)

onegene="NM_030930"
plotCounts(DEdds, gene=onegene, intgroup=c( "group"))