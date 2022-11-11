#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Rsubread")
library("Rsubread")

#set your out directory and the gtf for your genome
outdir="/Shares/down/heatshock/analysis_June2022/outfiles/PRO/gene_counts_and_diffexp_removechr/body_only/"
annot_file <-"/scratch/Shares/dowell/genomes/hg38/PROparts_genenames/Nov2022_removechrnotinfasta/body.saf"


setwd(outdir)

filetablename="/Shares/down/heatshock/analysis_June2022/scripts/T21HSscripts/metadata/PRObams_meta.txt"
filetable <- read.csv(filetablename, header=TRUE, sep="\t")
filelist <-as.vector(filetable$bamfile)


#check if all the bam files exist
if (!all(file.exists(filelist))) {
  print("WARNING: Not all specified files exist")
}

coverage <- featureCounts(files=filelist,
                          annot.ext=annot_file,
                          isGTFAnnotationFile=FALSE,
                          isPairedEnd=FALSE,
                          nthreads=32, allowMultiOverlap=TRUE)

colnames(coverage$counts) <- filetable$label

time <- strsplit(as.character(Sys.time()), split = " ")[[1]][2]
time <- paste(strsplit(time, split = ":")[[1]], collapse = '')
time


fileroot<-paste0(outdir, "featureCounts_body_", time)

#you can save the whole session as a R image. I would not suggest using this file. 
save.image(paste0(fileroot, ".RData"))

write.csv(coverage$counts, paste(fileroot,".coverage.csv", sep=""))
write.csv(coverage$stat, paste(fileroot,".stat.csv", sep=""))
write.csv(coverage$annotation, paste(fileroot,".annotation.csv", sep=""))
write.csv(coverage$targets, paste(fileroot,".targets.csv", sep=""))



