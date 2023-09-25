#2023-09-22: this script was used for pan-cancer snATAC paper
#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)

require(tidyr)
require(dplyr)
set.seed(1234)

library(plyr)
library(dplyr)
library(tibble)
library(reshape)

library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

library(JASPAR2020)
library(TFBSTools)
library(patchwork)

atac=readRDS(paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/3.Merge_snATAC/Merge.vers.20230114/PanCan_object/',
'CellsSampled.v.20230117/Pancan_225Samples_600cellsPerCluster_200PCs.chromVar.20230117.rds.gz',sep=''))

x=Motifs(atac[['pancan_s']])
mot=as.data.frame(x@positions)
peaks=as.data.frame(rownames(atac@assays$pancan_s))
peaks$Peaks_souce="Merged_PanCan_snATAC"


colnames(peaks)[1]='Peak'
mot$motif_cooord=paste(mot$seqnames,mot$start,mot$end,sep="-")
in_peaks=StringToGRanges(peaks$Peak, sep = c("-", "-"))
out_peaks=StringToGRanges(mot$motif_cooord, sep = c("-", "-"))
olap=as.data.frame(findOverlaps(in_peaks,out_peaks))
#pairs=cbind(da_p[olap$queryHits,],olap$queryHits,peaks[olap$subjectHits,],olap$subjectHits)
pairs=cbind(peaks[olap$queryHits,],mot[olap$subjectHits,])
pairs_1=pairs %>% dplyr::select ('Peak','group_name','strand','score','motif_cooord')
colnames(pairs_1)[5]= 'motif_coord'
jaspar=read.table('/diskmnt/Projects/HTAN_analysis/snATAC/Signac/CCRCC/JASPAR2020_motifs.txt',sep='\t',
header=TRUE)

colnames(jaspar)[1]='group_name'
pairs_2=merge(pairs_1,jaspar,all.x=TRUE)
write.table(pairs_2,paste('out/Motifs_matched.PanCan_snATAC_merged.object.20230118.tsv',sep=''),sep='\t',
quote=FALSE,row.names=FALSE)
write.table(mot,paste('out/All_motifs_before_matching.PanCan_snATAC_merged.object.ObjCellsSampled.20230118.tsv',sep=''),sep='\t',
quote=FALSE,row.names=FALSE)