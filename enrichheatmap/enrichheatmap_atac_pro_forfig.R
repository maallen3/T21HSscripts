library(EnrichedHeatmap)
library(GenomicRanges)
library(tidyr)
library(ggplot2)
library(plyranges)
library(circlize)



#Set in and out directories

##set outdir
mainPROindir="/Shares/down/heatshock/analysis_June2022/outfiles/PRO/"
mainATACindir="/Shares/down/heatshock/analysis_June2022/outfiles/ATAC/"
metadataindir="/Shares/down/heatshock/analysis_June2022/scripts/T21HSscripts/metadata/"
mainoutDir="/Shares/down/heatshock/analysis_June2022/outfiles/mix/"
subDir = "PROseq_around_HSF1_Chipseq_motif/"
outdir=file.path(mainoutDir, subDir)
dir.create(outdir, showWarnings = FALSE)

##chip and motif files
motiffile = "/scratch/Shares/dowell/motifs/HOCOMOCO_HUMAN_v11_p1e-6_grch38/HSF1_HUMAN.H11MO.0.A.bed"
chipfile= "/Shares/down/heatshock/analysis_June2022/extra_indata/HSF1chipseq_lymphoblastoidBeds_101713_peaks.bed"
genefile= "/scratch/Shares/dowell/genomes/hg38/hg38_refseq_genenames.bed"
genetssfile = "/scratch/Shares/dowell/genomes/hg38/hg38_refseq_genenames_TSS_1000bp.bed"

#PRO and ATAC metadataATACmetadatafile=#read the metadata
metadataindir="/Shares/down/heatshock/analysis_June2022/scripts/T21HSscripts/metadata/"
ATACmetadatafile=paste(metadataindir,"ATACbams_meta.txt", sep="")
ATACmetadata <- read.csv(ATACmetadatafile, header=TRUE, sep="\t")
head(ATACmetadata)
PROmetadatafile=paste(metadataindir,"PRObams_meta.txt", sep="")
PROmetadata <- read.csv(PROmetadatafile, header=TRUE, sep="\t")
head(PROmetadata)

##mumerge bed file
ATACpeakfile <- paste(mainATACindir,"/peak_merge/ATACHMMRATACmerge_MUMERGE.bed",sep="")
PROpeakfile <- paste(mainPROindir,"/bidir_merge/Tfitall_mumerge_MUMERGE.bed",sep="")

#size factors from Deseq2
ATACsizefactors=paste(mainATACindir,"/peaks_counts_and_diff/sizefactors.csv", sep="")
sfdf = read.csv(ATACsizefactors)


##set some standard paramaters
usecentermotif=TRUE
upwindowaroundmotif=500
downwindowaroundmotif=500
heatmapbinsize=50
ATACrepnum = "R1"
PROrepnum = "R2"


##create a metadata file to go with these plot
metadatafile=paste(outdir,"metadata.txt", sep="")
line = "This is the metadata for this directory."
write(line,file=metadatafile)
line = paste("motif file=", motiffile, sep="")
write(line,file=metadatafile,append=TRUE)
line = paste("PROpeakfile file=", PROpeakfile, sep="")
write(line,file=metadatafile,append=TRUE)
line = paste("chipfile file=", chipfile, sep="")
write(line,file=metadatafile,append=TRUE)
line = paste("ATACpeakfile file=", ATACpeakfile, sep="")
write(line,file=metadatafile,append=TRUE)
line = paste("usecentermotif=", usecentermotif, sep="")
write(line,file=metadatafile,append=TRUE)
line = paste("upwindowaroundmu=", upwindowaroundmotif, sep="")
write(line,file=metadatafile,append=TRUE)
line = paste("downwindowaroundmu=", downwindowaroundmotif, sep="")
write(line,file=metadatafile,append=TRUE)
line = paste("heatmapbinsize=", heatmapbinsize, sep="")
write(line,file=metadatafile,append=TRUE)
line = paste("ATACrepnum=", ATACrepnum, sep="")
write(line,file=metadatafile,append=TRUE)
line = paste("PROrepnum=", PROrepnum, sep="")
write(line,file=metadatafile,append=TRUE)

#import the motifs
motifs = read.table(motiffile)
motifs["loc"] = paste(motifs$V1,":",  motifs$V2,"-",motifs$V3, sep="")
grmotifs <- GRanges(seqnames = motifs[[1]], ranges = IRanges(motifs[[2]], motifs[[3]]), score = motifs[[5]], name= motifs[[4]],  strand=motifs[[6]], loc=motifs[,"loc"])
#when testing only look at chr1 to make this script faster
#grmotifs <- grmotifs %>% filter(seqnames=="chr1")

genes = read.table(genefile)
genes["loc"] = paste(genes$V1,":",  genes$V2,"-",genes$V3, sep="")
grgenes <- GRanges(seqnames = genes[[1]], ranges = IRanges(genes[[2]], genes[[3]]), score = genes[[5]], name= genes[[4]],  strand=genes[[6]], loc=genes[,"loc"])
#when testing only look at chr1 to make this script faster
#grmotifs <- grmotifs %>% filter(seqnames=="chr1")

genestss = read.table(genetssfile)
genestss["loc"] = paste(genestss$V1,":",  genestss$V2,"-",genestss$V3, sep="")
grgenestss <- GRanges(seqnames = genestss[[1]], ranges = IRanges(genestss[[2]], genestss[[3]]), score = genestss[[5]], name= genestss[[4]],  strand=genestss[[6]], loc=genestss[,"loc"])

#grab the center of the motifs and check if motifs exist on both strands in the same position
if(usecentermotif==TRUE){
  grmotifs=  mutate(anchor_center(grmotifs), width = 1)}
n_occur <- data.frame(table(motifs$loc))
bothstrands <- motifs[motifs$loc %in% n_occur$Var1[n_occur$Freq > 1],]
onestrand <- motifs[motifs$loc %in% n_occur$Var1[n_occur$Freq == 1],]
grmotifsbothstrands <- GRanges(seqnames = bothstrands[[1]], ranges = IRanges(bothstrands[[2]], bothstrands[[3]]), score = bothstrands[[5]], name= bothstrands[[4]],  strand=bothstrands[[6]], loc=bothstrands[,"loc"])
grmotifsonestrand <- GRanges(seqnames = onestrand[[1]], ranges = IRanges(onestrand[[2]], onestrand[[3]]), score = onestrand[[5]], name=onestrand[[4]],  strand=onestrand[[6]], loc=onestrand[,"loc"])


#write what we learned to the metadata file
line = paste("motif file, n=", length(grmotifs))
write(line,file=metadatafile,append=TRUE)
line = paste("There are , n=", length(grmotifsbothstrands)/2, " motifs on both strands. Therefore they are n=",length(grmotifsbothstrands), "of the total motifs.")
write(line,file=metadatafile,append=TRUE)


#I'm going to require my motifs to be in a chipseq region and a ATAC-seq peak

peaks = read.table(chipfile)
peaks2 = read.table(ATACpeakfile)
grpeakCHIP <- GRanges(seqnames = peaks[[1]], ranges = IRanges(peaks[[2]], peaks[[3]]))
grpeakATAC <- GRanges(seqnames = peaks2[[1]], ranges = IRanges(peaks2[[2]], peaks2[[3]]))
motifs_in_peak = join_overlap_inner(grmotifs, grpeakCHIP) # includes double the window above
motifs_in_peak = granges(motifs_in_peak)
motifs_in_peak = mutate(anchor_center(motifs_in_peak), width = 1)
motifs_in_peak <- GRanges(seqnames = seqnames(motifs_in_peak), ranges = IRanges(start(motifs_in_peak), start(motifs_in_peak)))
motifs_in_peak

motifs_in_peak_not_in_gene <- subsetByOverlaps(motifs_in_peak, grgenes, invert=TRUE)
motifs_in_peak_in_gene <- subsetByOverlaps(motifs_in_peak, grgenes)
motifs_in_peak<-motifs_in_peak_not_in_gene


#motifs_in_peak<-motifs_in_peak_in_gene

#write what we learned to the metadata file
line = paste("chip-seq or atac-seq regions, n=", length(grpeak))
write(line,file=metadatafile,append=TRUE)
line = paste("motifs in peaks, n=", length(motifs_in_peak))
write(line,file=metadatafile,append=TRUE)


#create matrix objects for each atac sample
matlist <- list()
quants <- list()

onerepATAC = ATACmetadata[ATACmetadata$replicate==ATACrepnum,]
onerepPRO = PROmetadata[PROmetadata$replicate==PROrepnum,]
atacsamplenames = as.vector(onerepATAC$group)
atacsamples = as.vector(onerepATAC$bedgraphfile)
atacsamplelabel = as.vector(onerepATAC$label)
prosamplenames = as.vector(onerepPRO$group)
prosamples = as.vector(onerepPRO$bedgraphfile)
prosamplelabel = as.vector(onerepPRO$label)


matposlist <- list()
matneglist <- list()
proquants <- list()


for (i in 1:length(prosamplenames)){
  message(i)
  samplename = prosamplenames[i]
  proseqfile<-prosamples[i]
  #read the bedgraph and split to postive and negative strand counts
  df = read.table(proseqfile)
  gr = GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]), coverage = df[[4]])
  pos = gr %>% filter(coverage>0)
  neg = gr %>% filter(coverage<0)
  mcols(neg)$coverage <-mcols(neg)$coverage*-1
  numberofregions = as.character(length(motifs_in_peak))
  #pos <= pos %>% filter(seqnames=="chr21")
  mat0pos <- normalizeToMatrix(pos, motifs_in_peak, value_column = "coverage",extend = c(upwindowaroundmotif,downwindowaroundmotif), w=heatmapbinsize)
  matposlist [[i]]<-mat0pos
  #neg <= neg %>% filter(seqnames=="chr21")
  mat0neg <- normalizeToMatrix(neg, motifs_in_peak, value_column = "coverage", extend = c(downwindowaroundmotif, upwindowaroundmotif),  w=heatmapbinsize)
  matneglist [[i]]<-mat0neg
  theseproquants = c(quantile(mat0pos, c(0.99)), quantile(mat0neg, c(0.99)))
  max99 <- max(unlist(theseproquants))
  proquants [i]<-max99
  #put info in the metadata file
  line = paste("bedgraph file=", proseqfile, sep="")
  write(line,file=metadatafile,append=TRUE)
  line = paste("proseq bedgraph file, lines=", length(gr))
  write(line,file=metadatafile,append=TRUE)
}




quantile(matposlist[[1]])
quantile(matposlist[[2]])
quantile(matposlist[[3]])
quantile(matposlist[[4]])

quantile(matneglist[[1]])
quantile(matneglist[[2]])
quantile(matneglist[[3]])
quantile(matneglist[[4]])

matlist <- list()
quants <- list()



for (i in 1:length(atacsamples)){
  #collect information for one SAMPLE
  message(i)
  bedgraphfilename=atacsamples[i]
  samplename = atacsamplenames[i]
  print(samplename)
  samplelabel = atacsamplelabel[i]
  normfactor = as.numeric(sfdf[sfdf$X==samplelabel,]["sf"])
  print(normfactor)
  #read the bedgraph AND normalize via Deseq2 size factors
  df = read.table(bedgraphfilename)
  df[[4]] <- df[[4]]/normfactor
  gr = GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]), coverage = df[[4]])
  #when testing only look at chr1 to make this script faster
  #gr <- gr %>% filter(seqnames=="chr1")
  numberofregions = as.character(length(motifs_in_peak))
  mat0pos <- normalizeToMatrix(gr, motifs_in_peak, value_column = "coverage",extend = c(upwindowaroundmotif,downwindowaroundmotif), w=heatmapbinsize)
  matlist[[i]]<-mat0pos
  max99 <- quantile(mat0pos, c(0.99))
  quants[i]<-max99
  #write what I did to the metadataout file
  line = paste("bedgraph file=", bedgraphfilename, sep="")
  write(line,file=metadatafile,append=TRUE)
  line = paste("ATACseq bedgraph file, lines=", length(gr))
  write(line,file=metadatafile,append=TRUE)
}

#______ graphing section



quantile(matlist[[1]])
quantile(matlist[[2]])
quantile(matlist[[3]])
quantile(matlist[[4]])


maxquant <- max(unlist(proquants))
maxquantatac <- max(unlist(quants))
#when max quant is 0 it has to be set manualy
maxquant<-.8
col_fun_fwd = colorRamp2(c(0, maxquant), c("white", "blue"))
col_fun_rev = colorRamp2(c(0, maxquant), c("white", "red"))
col_fun_atac = colorRamp2(c(0, maxquantatac), c("white", "brown"))

prominabsolute=0
promaxabsolute=0.5


colorlist=c('#32CD32', '#006400', '#DA70D6','#4B0082')
#run the first sample and make a list
samplename = prosamplenames[1]
mat0pos <- matposlist[[1]]
mat0neg <- matneglist[[1]]
if (promaxabsolute==0){
  p<-EnrichedHeatmap(mat0pos, col = col_fun_fwd, name=paste(samplename, "PRO+", sep=""), column_title = paste(samplename," +"))
  n<-EnrichedHeatmap(mat0neg, col = col_fun_rev, name=paste(samplename, "PRO-", sep=""), column_title = paste(samplename," -"))
}else{
  p<-EnrichedHeatmap(mat0pos, col = col_fun_fwd, name=paste(samplename, "PRO+", sep=""), column_title = paste(samplename," +"), top_annotation = HeatmapAnnotation(enriched = anno_enriched(ylim = c(prominabsolute, promaxabsolute),gp = gpar(col =colorlist[1]))))
  n<-EnrichedHeatmap(mat0neg, col = col_fun_rev, name=paste(samplename, "PRO-", sep=""), column_title = paste(samplename," -"), top_annotation = HeatmapAnnotation(enriched = anno_enriched(ylim = c(prominabsolute, promaxabsolute),gp = gpar(col =colorlist[1]))))
}

ht_list_pos <- p
ht_list_neg <- n




#loop through the other samples
for (i in 2:length(prosamplenames)){
  samplename = prosamplenames[i]
  message(i)
  mat0pos <- matposlist[[i]]
  mat0neg <- matneglist[[i]]
  if (promaxabsolute==0){
    p<-EnrichedHeatmap(mat0pos, col = col_fun_fwd, name=paste(samplename, "PRO+", sep=""), column_title = paste(samplename," +"))
    n<-EnrichedHeatmap(mat0neg, col = col_fun_rev, name=paste(samplename, "PRO-", sep=""), column_title = paste(samplename," -"))
  }else{
    p<-EnrichedHeatmap(mat0pos, col = col_fun_fwd, name=paste(samplename, "PRO+", sep=""), column_title = paste(samplename," +"), top_annotation = HeatmapAnnotation(enriched = anno_enriched(ylim = c(prominabsolute, promaxabsolute),gp = gpar(col =colorlist[i]))),show_heatmap_legend = FALSE)
    n<-EnrichedHeatmap(mat0neg, col = col_fun_rev, name=paste(samplename, "PRO-", sep=""), column_title = paste(samplename," -"), top_annotation = HeatmapAnnotation(enriched = anno_enriched(ylim = c(prominabsolute, promaxabsolute),gp = gpar(col =colorlist[i]))),show_heatmap_legend = FALSE)}
  ht_list_pos <- ht_list_pos+p
  ht_list_neg <- ht_list_neg+n}

ht_list <- ht_list_neg+ht_list_pos

draw(ht_list, merge_legends=TRUE, heatmap_legend_side="bottom", 
     annotation_legend_side="right")


D21diffposmatrix <- matposlist[[2]] - matposlist[[1]]
T21diffposmatrix <- matposlist[[4]] - matposlist[[3]]
D21diffnegmatrix <- matneglist[[2]] - matneglist[[1]]
T21diffnegmatrix <- matneglist[[4]] - matneglist[[3]]


col_fun_pro_diff <-colorRamp2(c(-2, 0, 2), c("orange", "white", "black"))


mindiff_pro=0
#if you don't know what to set maxdiff at use 0 and it will autoset but it will not be the same for all samples diffs
maxdiff_pro=0.5

if (maxdiff_pro==0){
  T21diffpos <- EnrichedHeatmap(T21diffposmatrix, name="diff T21 pos", col = col_fun_pro_diff,  column_title = "diff T21 pos")
  D21diffpos <- EnrichedHeatmap(D21diffposmatrix, name="diff D21 pos", col = col_fun_pro_diff, column_title = "diff D21 pos")
  T21diffneg <- EnrichedHeatmap(T21diffnegmatrix, name="diff T21 neg", col = col_fun_pro_diff,  column_title = "diff T21 neg")
  D21diffneg <- EnrichedHeatmap(D21diffnegmatrix, name="diff D21 neg", col = col_fun_pro_diff, column_title = "diff D21 neg")
}else{
  T21diffpos <- EnrichedHeatmap(T21diffposmatrix, name="diff T21 pos", col = col_fun_pro_diff,  column_title = "diff T21 pos", top_annotation =  HeatmapAnnotation(enriched = anno_enriched(ylim = c(mindiff_pro, maxdiff_pro),gp = gpar(col =colorlist[4]))))
  D21diffpos <- EnrichedHeatmap(D21diffposmatrix, name="diff D21 pos", col = col_fun_pro_diff, column_title = "diff D21 pos", top_annotation =  HeatmapAnnotation(enriched = anno_enriched(ylim = c(mindiff_pro, maxdiff_pro),gp = gpar(col =colorlist[2]))),show_heatmap_legend = FALSE)
  T21diffneg <- EnrichedHeatmap(T21diffnegmatrix, name="diff T21 neg", col = col_fun_pro_diff,  column_title = "diff T21 neg", top_annotation =  HeatmapAnnotation(enriched = anno_enriched(ylim = c(mindiff_pro, maxdiff_pro),gp = gpar(col =colorlist[4]))),show_heatmap_legend = FALSE)
  D21diffneg <- EnrichedHeatmap(D21diffnegmatrix, name="diff D21 neg", col = col_fun_pro_diff, column_title = "diff D21 neg", top_annotation =  HeatmapAnnotation(enriched = anno_enriched(ylim = c(mindiff_pro, maxdiff_pro),gp = gpar(col =colorlist[2]))),show_heatmap_legend = FALSE)}


ht_list_wtih_diffs = D21diffneg+D21diffpos+T21diffneg+T21diffpos
draw(ht_list_wtih_diffs, main_heatmap="diff D21 pos")

atacminabsolute=0
atacmaxabsolute=16

samplename = atacsamplenames[1]
mat <- matlist[[1]]
if (atacmaxabsolute==0){
  b<-EnrichedHeatmap(mat, col = col_fun_atac, name=paste(samplename, " ATAC", sep=""), column_title = paste(samplename," +"))
}else{
  b<-EnrichedHeatmap(mat, col = col_fun_atac, name=paste(samplename, "ATAC", sep=""), column_title = paste(samplename," +"), top_annotation = HeatmapAnnotation(enriched = anno_enriched(ylim = c(atacminabsolute, atacmaxabsolute),gp = gpar(col =colorlist[1]))))}
ht_list_atac <- b

#loop through the other samples
for (i in 2:length(atacsamplenames)){
  samplename = atacsamplenames[i]
  message(i)
  mat <- matlist[[i]]
  if (atacmaxabsolute==0){
    b<-EnrichedHeatmap(mat, col = col_fun_atac, name=paste(samplename, " ATAC", sep=""), column_title = paste(samplename," +"))
  }else{
    b<-EnrichedHeatmap(mat, col = col_fun_atac, name=paste(samplename, " ATAC", sep=""), column_title = paste(samplename," +"), top_annotation = HeatmapAnnotation(enriched = anno_enriched(ylim = c(atacminabsolute, atacmaxabsolute),gp = gpar(col =colorlist[i]))),show_heatmap_legend = FALSE)}
  ht_list_atac <- ht_list_atac+b}

ht_list_atac


D21diffmatrix <- matlist[[2]] - matlist[[1]]
T21diffmatrix <- matlist[[4]] - matlist[[3]]

col_fun_atac_diff <-colorRamp2(c(-60, 0, 60), c("orange", "white", "black"))

mindiff_atac=0
#if you don't know what to set maxdiff at use 0 and it will autoset but it will not be the same for all samples diffs
maxdiff_atac=4

if (maxdiff_atac==0){
  T21diff <- EnrichedHeatmap(T21diffmatrix, name="diff T21", col = col_fun_atac_diff,  column_title = "diff T21")
  D21diff <- EnrichedHeatmap(D21diffmatrix, name="diff D21", col = col_fun_atac_diff, column_title = "diff D21")
  }else{
  T21diff <- EnrichedHeatmap(T21diffmatrix, name="diff T21", col = col_fun_atac_diff,  column_title = "diff T21", top_annotation = HeatmapAnnotation(enriched = anno_enriched(ylim = c(mindiff_atac, maxdiff_atac),gp = gpar(col =colorlist[4]))))
  D21diff <- EnrichedHeatmap(D21diffmatrix, name="diff D21", col = col_fun_atac_diff, column_title = "diff D21", top_annotation = HeatmapAnnotation(enriched = anno_enriched(ylim = c(mindiff_atac, maxdiff_atac),gp = gpar(col =colorlist[2]))),show_heatmap_legend = FALSE)
  }


#the most diff open chromatin in T21 orders it
htdiff <- T21diff+D21diff+ht_list_wtih_diffs
htdata <- ht_list_atac+ht_list
  
ht <- htdiff+ htdata

draw(ht, merge_legend=TRUE,main_heatmap="diff T21 pos")

ht2 <- D21diff+T21diff+ht_list_atac+ht_list_wtih_diffs+ht_list
draw(ht2, merge_legend=TRUE,main_heatmap="diff T21 pos",heatmap_legend_side="bottom", 
     annotation_legend_side="right")

svg(paste(outdir,"EnrichedHeatmap_of_ATAC_PRO_overHSFwithchip_nogenes.svg", sep=""))
draw(ht2, merge_legend=TRUE,main_heatmap="D21_fourtytwoPRO+",heatmap_legend_side="bottom", 
     annotation_legend_side="right")
dev.off()

draw(ht2, merge_legend=TRUE,main_heatmap="D21_fourtytwoPRO-",heatmap_legend_side="bottom", 
     annotation_legend_side="right")

draw(ht2, merge_legend=TRUE,main_heatmap="D21_fourtytwo ATAC",heatmap_legend_side="bottom", 
     annotation_legend_side="right")

#plot the row means for all positions in the matrix
results=lapply(matposlist, FUN=colMeans)
colmeansdf=as.data.frame(do.call(rbind, results))
colmeansdf=t(colmeansdf)
colnames(colmeansdf) <-prosamplenames
colmeansdf=as.data.frame(colmeansdf)
colmeansdf$position <- seq(-1*upwindowaroundmotif,downwindowaroundmotif-1,heatmapbinsize)
meandflong <- tidyr::gather(colmeansdf, condition, depth, -position) %>% tidyr::separate(condition, c("genotype", "tempature"))
meandflong$strand <- "+"

#plot the metaplot for PRO
#plot the row means for all positions in the matrix pos
results=lapply(matposlist, FUN=colMeans)
colmeansdf=as.data.frame(do.call(rbind, results))
colmeansdf=t(colmeansdf)
colnames(colmeansdf) <-prosamplenames
colmeansdf=as.data.frame(colmeansdf)
colmeansdf$position <- seq(-1*upwindowaroundmotif,downwindowaroundmotif-1,heatmapbinsize)
meandflong <- tidyr::gather(colmeansdf, condition, depth, -position) %>% tidyr::separate(condition, c("genotype", "tempature"))
meandflong$strand <- "+"

#plot the row means for all positions in the matrix neg
results1=lapply(matneglist, FUN=colMeans)
colmeansdf1=as.data.frame(do.call(rbind, results1))
colmeansdf1=t(colmeansdf1)
colnames(colmeansdf1) <-prosamplenames
colmeansdf1=as.data.frame(colmeansdf1)
colmeansdf1$position <- seq(-1*downwindowaroundmotif,upwindowaroundmotif-1,heatmapbinsize)
meandflong1 <- tidyr::gather(colmeansdf1, condition, depth, -position) %>% tidyr::separate(condition, c("genotype", "tempature"))
meandflong1$depth <- meandflong1$depth*-1
meandflong1$strand <- "-"

meandflong_bothstrands <- rbind(meandflong, meandflong1)
meandflong_bothstrands$label <- paste(meandflong_bothstrands$genotype,meandflong_bothstrands$tempature)

ggplot(data=meandflong_bothstrands, aes(x=position, y=depth,colour=label, linetype=interaction(genotype, tempature, strand))) +
  geom_line(size=2)+scale_linetype_manual(values=c("solid","solid", "dashed","dashed", "solid","solid", "dashed","dashed"))+scale_color_manual(values=c('#32CD32', '#DA70D6','#006400', '#4B0082', '#32CD32', '#DA70D6','#006400', '#4B0082'))+theme(text = element_text(size=20))+ theme_minimal()+ theme(text = element_text(size = 20))


svg(paste(outdir,"Rowmeans_of_PRO_overHSFwithchip_no_genes.svg", sep=""))
ggplot(data=meandflong_bothstrands, aes(x=position, y=depth,colour=label,linetype=interaction(strand))) +geom_line(size=2)+scale_color_manual(values=c("D21 thirtyseven"='#32CD32', "D21 fourtytwo"='#006400',"T21 thirtyseven"='#DA70D6', "T21 fourtytwo"='#4B0082'))+scale_linetype_manual(values=c("+"="solid", "-"="dashed"))+theme(text = element_text(size=20))+ theme_minimal()+ theme(text = element_text(size = 20))
dev.off()

rowmeansfig <-ggplot(data=meandflong_bothstrands, aes(x=position, y=depth,colour=label,linetype=interaction(strand))) +geom_line(size=2)+scale_color_manual(values=c("D21 thirtyseven"='#32CD32', "D21 fourtytwo"='#006400',"T21 thirtyseven"='#DA70D6', "T21 fourtytwo"='#4B0082'))+scale_linetype_manual(values=c("+"="solid", "-"="dashed"))+theme(text = element_text(size=20))+ theme_minimal()+ theme(text = element_text(size = 20))
rowmeansfig
#ggsave(file=paste(outdir,"Rowmeans_of_PRO_overHSFwithchip.svg", sep=""), plot=rowmeansfig)


  