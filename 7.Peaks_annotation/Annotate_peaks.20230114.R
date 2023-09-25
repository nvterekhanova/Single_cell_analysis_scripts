#2023-09-22: this script was used for pan-cancer snATAC paper

#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp
library(Signac)
library(Seurat)
library(presto)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(reshape)
library(reshape2)
library(EnsDb.Hsapiens.v86)

peaks=read.table(paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/3.Merge_snATAC/Merge.vers.20230114/',
'Select_peaks.20230114/PanCan/peaks/PanCanrecentered_final.filtered.20230114.tsv',sep=''),sep='\t',header=T)

peaks_1=StringToGRanges(peaks$new_peak, sep = c("-", "-"))
###Now annotate peaks:
peakAnno <- annotatePeak(peaks_1, tssRegion=c(-1000, 100),
                         TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db")
                         
anno=as.data.frame(peakAnno)
peaks$Gene=anno$SYMBOL
peaks$Type=anno$annotation
peaks$geneId=anno$geneId
peaks$peak_distanceToTSS=anno$distanceToTSS
peaks$geneID=anno$geneID

dir.create('out')
write.table(peaks, "out/PanCanrecentered_final.filtered.Annotated.20230114.tsv",sep='\t',
row.names=F,quote=F)
