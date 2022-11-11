library("Rsubread")

metadatafile="/Shares/down/heatshock/analysis_June2022/scripts/T21HSscripts/metadata/PRObams_meta.txt"
outdir="/Shares/down/heatshock/analysis_June2022/outfiles/PRO/bidir_merge_diff/"
indir="/Shares/down/heatshock/analysis_June2022/outfiles/PRO/bidir_merge_diff/"
annot_file = paste(indir, "PRO_seq_MUMERGE.saf", sep="")

safdf = read.table(annot_file, sep=",", header=TRUE)
head(safdf)
tail(safdf)
dim(safdf)


filetable <- read.csv(metadatafile, header=TRUE, sep="\t")
#I have a column in my metadatafiles called "bamfile" that is the bam files
#I also have a column called "label"
filelist <-as.vector(filetable$bamfile)

#check if all the bam files exist
if (!all(file.exists(filelist))) {
  print("WARNING: Not all specified files exist")
}

coverage <- featureCounts(files=filelist,
                          annot.ext=safdf,
                          isGTFAnnotationFile=FALSE,
                          isPairedEnd=FALSE,
                          nthreads=32)

colnames(coverage$counts) <- filetable$label

time <- strsplit(as.character(Sys.time()), split = " ")[[1]][2]
time <- paste(strsplit(time, split = ":")[[1]], collapse = '')
time

#you can save the whole session as a R image. I would not suggest it. 
#save.image(paste0(outdir, "outfilename", GTFattrType, "full", "_", time, ".RData"))


fileroot<-paste0(outdir, "PRO_tfit_mumerge_peaks","_", time)
#fileroot<-paste0(outdir, "featureCounts_tss_", time)

save.image(paste0(fileroot, ".RData"))

write.csv(coverage$counts, paste(fileroot,".coverage.csv", sep=""))
write.csv(coverage$stat, paste(fileroot,".stat.csv", sep=""))
write.csv(coverage$annotation, paste(fileroot,".annotation.csv", sep=""))
write.csv(coverage$targets, paste(fileroot,".targets.csv", sep=""))
