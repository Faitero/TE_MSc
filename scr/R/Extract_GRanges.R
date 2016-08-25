### Extract data from GRanges
##################################
### overlap Deseq Log2fold with actual clasification  
##################################
setwd('/home/igor/Dropbox (UCL-MN Team)/hnRNPC.NMD (1).manuscript/Alus_primates_evolution/Data')  ## Curro
setwd('/Users/Igor/Dropbox (UCL-MN Team)/hnRNPC.NMD (1).manuscript/Alus_primates_evolution/Data') ## CAsa

library(GenomicRanges)
library(rtracklayer)

whole_final <-  read.csv('whole_final.tab', header = TRUE, sep = "\t")

whole_final_bed <- data.frame(whole_final$chr_hg19_t, whole_final$start_hg19_t, whole_final$end_hg19_t, whole_final$aluexon, whole_final$Exon_Originated, whole_final$strand_hg19_t)
write.table(whole_final_bed, file="whole_final.bed", sep="\t", col.names =FALSE, row.names=FALSE)

## Aluinformation.plotted.inFig5.gr.RData’ -> this is DEXseq only Aluexons but including random intronic Alus. Maybe better you work with this:(‘supplfile10_DEXseq.summary.RData’ -> this is all exons. I wanted to add this as supply file to the paper, so good to have for the Aluexons in here the AluexonID from cross-species comparison

load('supplfile10_DEXseq.summary.RData')

head(results.gr)
seqnames(results.gr)
seqinfo(results.gr)
results.gr


elementMetadata(results.gr)

df_DESEQ <- as.data.frame(results.gr)
df_DESEQ_IDS <- data.frame(df_DESEQ$seqnames, df_DESEQ$start, df_DESEQ$end, rownames(df_DESEQ), df_DESEQ$exon_ID, df_DESEQ$strand, df_DESEQ$exon.class.redone)
write.table(df_DESEQ_IDS, file="DF_DESEQ.bed", sep="\t", col.names =FALSE, row.names=FALSE)



##################################
### Extract data from Granges
################################## bedtools intersect -s -wb -a DF_DESEQ_sort_alu.exon_NOQUOTES.bed -b whole_final_sort_flank_NOQUOTES.bed > delete.bed

DF_DESEQ_aluID <- read.table('DF_SESEQ_alu.exon_aluID.bed', header = FALSE, sep = "\t")

names(DF_DESEQ_aluID) <- c("chr", "start", "end", "rowname_ID", "exon_ID", "strand", "chr3ss", "start3ss", "end3ss", "aluexon_ID", "Exon_Originated", "strand3ss")
df_DESEQ_complete <- merge.data.frame(df_DESEQ, DF_DESEQ_aluID, by="exon_ID", all.x=TRUE) 


df_DESEQ_complete <- df_DESEQ_complete[-c(40:47,50)]

write.table(df_DESEQ_complete, file='supplfile10_DEXseq.summary.withAluID_TEST.tab', sep="\t", col.names =TRUE, row.names=TRUE)


