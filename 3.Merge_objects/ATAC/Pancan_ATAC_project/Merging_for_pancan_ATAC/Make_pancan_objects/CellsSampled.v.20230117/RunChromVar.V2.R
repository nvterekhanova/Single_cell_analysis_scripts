####Important!!! --to limit number of cores:
library(BiocParallel)
register(SerialParam())
register(MulticoreParam(15)) #number of cores to use - otherwise it crushes

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)

require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)
set.seed(1234)

library(plyr)
library(dplyr)
library(tibble)
#library(reshape)
library(plyr)

library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

library(ggplot2)
library(RColorBrewer)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)

#atac=readRDS("Pancan_600cellsPerCluster_200PCs.20230117.rds.gz")
atac=combined
DefaultAssay(atac)='pancan_s'

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# add motif information
atac <- AddMotifs(
  object = atac,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

atac <- RunChromVAR(
  object = atac,
  genome = BSgenome.Hsapiens.UCSC.hg38
)


saveRDS(atac,"Pancan_225Samples_600cellsPerCluster_200PCs.chromVar.20230117.rds.gz",compress=T)


#atac <- RunUMAP(object = atac, reduction = 'lsi', dims = 2:75,reduction.name = "umap.75PC")



meta=read.table('All_137_samples_metadata_data_freeze_v2.0.tsv',sep='\t',header=T)
meta$V1=gsub('HT137B1_S17','HT137B1_S1H7',meta$V1)
meta$Piece_ID=gsub('HT137B1_S17','HT137B1_S1H7',meta$Piece_ID)
meta$cell_type_1=meta$cell_type.harmonized
meta$cell_type_1=ifelse(meta$cell_type.harmonized %in% c('Hepatocytes','Enterocytes','Ciliated','Alveolar','Goblet',
'Cholangiocytes','Islets','Acinar','Adipocytes','Erythrocytes','Lymphocytes','OPC','Unknown'),'Other',
meta$cell_type.harmonized)

cell_t_cols=c(brewer.pal(name='Dark2',n=8),brewer.pal(name='Set2',n=8),'#C0C0C0')
names(cell_t_cols)=unique(meta$cell_type_1)[c(1,3,4,2,5:10,17,12:16,11)]

orig_1=as.data.frame(atac$Piece_ID)
meta$Barcode=gsub('.*_.*_(.*)','\\1',meta$V1)
meta$Barcode=paste(meta$Piece_ID,meta$Barcode,sep='_')
meta$Barcode=paste(meta$Cancer,meta$Barcode,sep='_')
rownames(meta)=meta$Barcode

meta=meta[rownames(orig_1),]

atac$cell_type.harmonized=meta$cell_type.harmonized
atac$cell_type_1=meta$cell_type_1

#####Subset Tumor cells:
ATAC=subset(atac,cell_type.harmonized=='Tumor')




#########ENDS HERE FOR NOW##########
####Later part for saving motifs:###
####################################

atac=readRDS('26_ccRCC_snATAC.selectedPeaks.chromvar.v3.20210602.rds')

x=Motifs(atac[['peaksMACS2']])
mot=as.data.frame(x@positions)
peaks=as.data.frame(rownames(atac@assays$peaksMACS2))
peaks$Peaks_souce="Merged_snATAC"


colnames(peaks)[1]='Peak'
mot$motif_cooord=paste(mot$seqnames,mot$start,mot$end,sep="-")
in_peaks=StringToGRanges(peaks$Peak, sep = c("-", "-"))
out_peaks=StringToGRanges(mot$motif_cooord, sep = c("-", "-"))
olap=as.data.frame(findOverlaps(in_peaks,out_peaks))
#pairs=cbind(da_p[olap$queryHits,],olap$queryHits,peaks[olap$subjectHits,],olap$subjectHits)
pairs=cbind(peaks[olap$queryHits,],mot[olap$subjectHits,])
pairs_1=pairs %>% dplyr::select ('Peak','group_name','strand','score','motif_cooord')
colnames(pairs_1)[5]= 'motif_coord'
jaspar=read.table('/diskmnt/Projects/HTAN_analysis/snATAC/Signac/CCRCC/JASPAR2020_motifs.txt',sep='\t',header=TRUE)
colnames(jaspar)[1]='group_name'
pairs_2=merge(pairs_1,jaspar,all.x=TRUE)
write.table(pairs_2,paste('Motifs_matched.26_snATAC_merged.object.20210531.tsv',sep=''),sep='\t',
quote=FALSE,row.names=FALSE)
write.table(mot,paste('All_motifs_before_matching.26_snATAC_merged.object.20210531.tsv',sep=''),sep='\t',
quote=FALSE,row.names=FALSE)
peaks_sel=read_delim('../../6.DA_motifs/Match_motifs/20210520/DEG_associated_Peaks.20210517.v1.tsv',delim='\t')
peaks_sel=as.data.frame(peaks_sel)
peaks_sel_1=merge(peaks_sel,pairs_2,all.x=TRUE)
write.table(peaks_sel_1, 'Motifs_matched.DEG_associated_Peaks.20210517.v1.tsv',sep="\t",row.names=FALSE,quote=FALSE)
