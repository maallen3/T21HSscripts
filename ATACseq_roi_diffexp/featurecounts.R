#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Rsubread")
library("Rsubread")


#turn mumerge file into saf
mumergedir="/Shares/down/heatshock/analysis_June2022/outfiles/ATAC/peak_merge/"
mumergefile=paste(mumergedir, "ATACHMMRATACmerge_MUMERGE.bed", sep="")
df = read.table(mumergefile, sep="\t", col.name=c("Chr", "Start", "End"))
#add a geneID and strand column. In this case strand doesn't matter so I use a .
df$GeneID = paste(df$Chr, ":", df$Start, "-", df$End, sep="")
df$Strand = "+"
#reanage the columns to standard saf format
df <- df %>% dplyr::select(GeneID, Chr, Start, End, Strand)
#output the saf file
write.table(df, paste(mumergedir, "ATACHMMRATACmerge_MUMERGE.saf", sep=""), quote = FALSE, row.names = FALSE, sep="\t")

#set your out directory and the gtf for your genome
outdir="/Shares/down/heatshock/analysis_June2022/outfiles/ATAC/peaks_counts_and_diff/"
annot_file <-paste(mumergedir, "ATACHMMRATACmerge_MUMERGE.saf", sep="")


setwd(outdir)

filetablename="/Shares/down/heatshock/analysis_June2022/scripts/T21HSscripts/metadata/ATACbams_meta.txt"
filetable <- read.csv(filetablename, header=TRUE, sep="\t")
filelist <-as.vector(filetable$bamfile)


#check if all the bam files exist
if (!all(file.exists(filelist))) {
  print("WARNING: Not all specified files exist")
}

coverage <- featureCounts(files=filelist,
                          isPairedEnd=TRUE,
                          nthreads=32,
                          annot.ext=annot_file)

colnames(coverage$counts) <- filetable$label

time <- strsplit(as.character(Sys.time()), split = " ")[[1]][2]
time <- paste(strsplit(time, split = ":")[[1]], collapse = '')
time


fileroot<-paste0(outdir, "featureCounts_open_", time)

#you can save the whole session as a R image. I would not suggest using this file. 
save.image(paste0(fileroot, ".RData"))

write.csv(coverage$counts, paste(fileroot,".coverage.csv", sep=""))
write.csv(coverage$stat, paste(fileroot,".stat.csv", sep=""))
write.csv(coverage$annotation, paste(fileroot,".annotation.csv", sep=""))
write.csv(coverage$targets, paste(fileroot,".targets.csv", sep=""))



