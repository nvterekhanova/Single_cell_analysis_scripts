export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp

library(Signac)
library(Seurat)
library(plyr)
library(dplyr)

dir_with_RDS_objs=''
path_to_ATAC_meta_data=''

cancers=c('HNSCC','CRC','BRCA','OV','CESC','GBM','ccRCC','PDAC','MM','UCEC')
for (disease in cancers){
print(disease)
atac=readRDS(paste(dir_with_RDS_objs,'/', disease,'_snATAC_Merged.PancanSet.','noDoublets.GeneActivity.20211024.rds',sep=''))

annot=read.table(path_to_ATAC_meta_data, sep='\t', header=T)
cell_types=unique(annot$cell_type.harmonized.cancer)

annot_sel=annot[annot$Cancer==disease,]
orig_1=as.data.frame(atac$Piece_ID)
rownames(annot_sel)=annot_sel$V1
annot_sel=annot_sel[rownames(orig_1),]
atac$cell_type.harmonized.cancer=annot_sel$cell_type.harmonized.cancer

atac$ID=paste(atac$cell_type.harmonized.cancer, atac$Piece_ID,sep='_')
all_ids=table(atac$ID)
sel_ids=names(all_ids[all_ids>50])
atac$ID=ifelse(atac$ID %in% sel_ids, atac$ID, 'Other')

ATAC=atac
DefaultAssay(ATAC)<-'ATACGeneActivity'
Idents(ATAC)=ATAC$ID

var_genes=row.names(ATAC)

aliquot.averages <- AverageExpression(ATAC, assays = 'ATACGeneActivity', slot ='data',features=var_genes)

file2write <- paste0("out/PieceID_cellType/",disease,"_AverageGeneActivity.byPieceID_cellType.AllGenes.",
"20220103.tsv")
write.table(aliquot.averages, file = file2write, quote = F, sep = "\t", row.names = T)
}

