
outdir="/Shares/down/heatshock/analysis_June2022/outfiles/PRO/bidir_merge_diff/"

mumergefile<-"/Shares/down/heatshock/analysis_June2022/outfiles/PRO/bidir_merge/Tfitall_mumerge_MUMERGE.bed"

df = read.table(mumergefile, sep="\t", col.name=c("Chr", "Start", "End"))

#add a geneID and strand column. In this case strand doesn't matter so I use a .
df$GeneID = paste(df$Chr, ":", df$Start, "-", df$End, sep="")
df$Strand = "."


#reanage the columns to standard saf format
df <- df %>% dplyr::select(GeneID, Chr, Start, End, Strand)

#output the saf file
write.csv(df, paste(outdir, "PRO_seq_MUMERGE.saf", sep=""), quote = FALSE, row.names = FALSE)