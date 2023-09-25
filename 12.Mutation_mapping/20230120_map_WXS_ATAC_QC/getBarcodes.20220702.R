#Based on getAssay.20220626.R
#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp
library(Signac)
library(Seurat)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(reshape)
library(reshape2)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(reshape)


disease='PDAC'
atac_tab=read.table(paste('ATAC_catalog.BAM_paths.20230120.v1.txt',sep=''),sep='\t',header=T)

atac_tab_s=atac_tab[atac_tab$Disease.Type==disease & atac_tab$Sample.Type %in% c('Tumor','Met','Normal'),]


for (piece_id in atac_tab_s$Piece_ID){
print(piece_id)
sample=atac_tab_s$Sample.ID[atac_tab_s$Piece_ID==piece_id]
data_type=atac_tab_s$Data.Type[atac_tab_s$Piece_ID==piece_id]
if (data_type=='snATAC'){
path=paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/0.Individual_objects/snATAC/out/',
disease,'/',disease,'_',sample,'/',disease,'_',sample,'_processed_atac.rds',sep='')
}else if (data_type=='10x_SC_Multi_ATAC_SEQ'){
path=paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/0.Individual_objects/Multiome/out/',
disease,'/',disease,'_',sample,'/',disease,'_',sample,'_processed_multiome.rds',sep='')
}
atac=readRDS(path)

barc=as.data.frame(colnames(atac))
write.table(barc,paste('out/',disease,'/barcodes/',piece_id,'.ATAC_barcodes.tsv',sep=''),sep='\t',quote=FALSE,
row.names=F,col.names=F)
}
