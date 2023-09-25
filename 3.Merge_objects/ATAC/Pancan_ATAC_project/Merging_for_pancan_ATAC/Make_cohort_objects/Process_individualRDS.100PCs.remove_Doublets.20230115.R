#2023-09-22: this script was used for pan-cancer snATAC paper
#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(RColorBrewer)

require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
set.seed(1234)

library(plyr)
library(dplyr)
library(tibble)

library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(ggplot2)
library(future)

###some parallelization-solution from the tutorial:
plan("multiprocess", workers = 20)
options(future.globals.maxSize = 200 * 1024^3) # for 200 Gb RAM


args = commandArgs(trailingOnly=TRUE)
disease=args[1]
date='20230114'
date_2='20230114'

con<-file(paste('logs/Remove_doublets/',disease,'_Remove_doublets.log.',date_2,'.txt',sep=''))
sink(con,append=TRUE)
sink(con,type='message',append=TRUE)

doublets_tab=read_delim('Doublet_annotation.v.20230114/PanCan/Doublet_Barcode_annotation.20230114.tsv',delim='\t')
doublets_tab=as.data.frame(doublets_tab)

atac=readRDS(paste('../Quantify_PanCanSet_perCohort/RDS/',disease,'_snATAC_Merged.PancanSet.',date,'.rds',sep=''))

###Change barcodes for PKD-sample (Need that for one PKD-sample); and also update barcodes for BRCA samples
if (disease=='ccRCC'){
   orig_1=as.data.frame(atac$Piece_ID)
   orig_1$Barcode=rownames(orig_1)
   orig_1$Barcode_new=gsub('^PKD_','ccRCC_',orig_1$Barcode)

   atac=RenameCells(atac,new.names=orig_1$Barcode_new)
   atac$Piece_ID=gsub('PKD_K1900070','ccRCC_K1900070',atac$Piece_ID)

#doublets_cancer$Disease=ifelse(doublets_cancer$Sample_ID=='TWAW-K1900070_1FB-XBa1','PKD',doublets_cancer$Disease)
}


doublets_cancer=doublets_tab[doublets_tab$Disease==disease,]
rownames(doublets_cancer)=paste(doublets_cancer$Disease,doublets_cancer$Piece_ID, doublets_cancer$Barcodes,sep='_')

orig_1=as.data.frame(atac$Piece_ID)
doublets_cancer=doublets_cancer[rownames(orig_1),]
atac$Doublet_annot=doublets_cancer$Doublet_final

#make a check:
#table(atac@meta.data[,c('Doublet_annot','Piece_ID')])

atac=subset(atac, Doublet_annot!='TRUE' | is.na(Doublet_annot))

DefaultAssay(atac)<-'pancan'

atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = 'q0')
atac <- RunSVD(atac,n=100)

atac <- RunUMAP(object = atac, reduction = 'lsi', dims = 2:50)
atac <- FindNeighbors(object = atac, reduction = 'lsi', dims = 2:50)


###this helps:
options(future.globals.maxSize= 891289600)
atac <- FindClusters(object = atac, verbose = FALSE, algorithm = 3)
atac <- FindClusters(object = atac, verbose = FALSE, algorithm = 3,resolution=2)

atac[['peaksinters']]<-NULL

meta.data=as.data.frame(atac@meta.data)
embeddings=as.data.frame(Embeddings(atac, reduction = "umap")[,1:2])
embeddings=embeddings[rownames(meta.data),]
meta.data$UMAP_1=embeddings$UMAP_1
meta.data$UMAP_2=embeddings$UMAP_2
write.table(meta.data, paste("RDS.withUMAPs.50PCs.noDoublets/Meta_data.",date_2,"/",disease,
"_snATAC_Merged.PancanSet.noDoublets.MetaData",date_2,".tsv",sep=''),sep='\t',quote=F,row.names=T)

saveRDS(atac, paste("RDS.withUMAPs.50PCs.noDoublets/",disease,"_snATAC_Merged.PancanSet.noDoublets.",date_2,
".rds",sep=''))


p1 <- DimPlot(atac, group.by = 'Piece_ID', pt.size = 0.1,label=T) +
ggplot2::ggtitle(paste(disease," snATAC samples",sep=''))
pdf(paste("RDS.withUMAPs.50PCs.noDoublets/plots/",disease,"PanCan_PeakSet_UMAP_50PCs.noDoublets.",date_2,".pdf",
sep=""),height=6,width=8)
p1
dev.off()

Idents(atac)=atac$pancan_snn_res.0.8
p1 <- DimPlot(atac,pt.size = 0.1,label=T) + ggplot2::ggtitle(paste(disease," snATAC samples",sep=''))
pdf(paste("RDS.withUMAPs.50PCs.noDoublets/plots/Seurat_clusters/",disease,
"PanCan_PeakSet_UMAP_50PCs.noDoublets.res.0.8.",date_2,".pdf",sep=""),height=6,width=8)
p1
dev.off()

Idents(atac)=atac$pancan_snn_res.2
p1 <- DimPlot(atac,pt.size = 0.1,label=T) + ggplot2::ggtitle(paste(disease," snATAC samples",sep=''))
pdf(paste("RDS.withUMAPs.50PCs.noDoublets/plots/Seurat_clusters/",disease,
"PanCan_PeakSet_UMAP_50PCs.noDoublets.res.2.",date_2,".pdf",sep=""),height=6,width=8)
p1
dev.off()


#####ENDS HERE
