#2023-09-22: script to generate average expression matrix per cell group. This script was used for pan-cancer snATAC project.
#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp
library(Signac)
library(Seurat)
library(sctransform)
library(glmGamPoi)

###use snRNA obj:

data_dir='/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v7.0'

rna=readRDS(paste(data_dir,'/snRNA/PanCan_merge/',
'Merged_206_snRNA.SCTrandformed.3KVarFeatures.TumorNormal.v7.20230118.rds',sep=''))

meta=read.table(paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/DATA/Objects.v7/',
'All_206_RNA_samples_metadata_data_freeze_v7.0.2.withDetailed.tsv',sep=''),sep='\t',header=T)
rownames(meta)=meta$X
meta$cell_type.detailed=ifelse(!is.na(meta$cell_type.normal),meta$cell_type.normal,
meta$cell_type.harmonized.cancer)

orig_1=as.data.frame(rna$Piece_ID)
meta=meta[rownames(orig_1),]
rna$cell_type.detailed=meta$cell_type.detailed



rna$Disease=ifelse(rna$Piece_ID %in% c("HT268B1-Th1H3", "HT029B1-S1PC", "HT035B1-S1PA",
 "HT1408-06","HT141B1-S1H1", "HT206B1-S1H4", "HT271B1-S1H3",
"HT378B1-S1H1", "HT378B1-S1H2", "HT384B1-S1H1", "HT517B1-S1H1"), "BRCA_Basal", rna$Disease)

rna$Piece_ID=paste(rna$Disease,rna$Piece_ID, sep='_')

rna$ID=paste(rna$Piece_ID,rna$cell_type.detailed,sep='__')

var_genes <- rownames(rna)
DefaultAssay(rna)<-'SCT'

Idents(rna)<-rna$ID
aliquot.averages <- AverageExpression(rna, assays = 'SCT', slot ='data',features=var_genes)

file2write <- paste0("out/AverageExpressionSCT_SlotData_AllGenes.20230125.tsv")
write.table(aliquot.averages, file = file2write, quote = F, sep = "\t", row.names = T)

n_cells=as.data.frame(table(rna$ID))
colnames(n_cells)=c('Cell_type','Cell_count')
write.table(n_cells,'out/N_cells_perGroup.RNA.20230125.tsv',sep='\t',quote = F, row.names =F)


#########################################
####Also do Normal cell type combined.###
#########################################

data_dir='/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v7.0'

rna=readRDS(paste(data_dir,'/snRNA/PanCan_merge/',
'Merged_206_snRNA.SCTrandformed.3KVarFeatures.TumorNormal.v7.20230118.rds',sep=''))

meta=read.table(paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/DATA/Objects.v7/',
'All_206_RNA_samples_metadata_data_freeze_v7.0.2.withDetailed.tsv',sep=''),sep='\t',header=T)
rownames(meta)=meta$X
meta$cell_type.detailed=ifelse(!is.na(meta$cell_type.normal),meta$cell_type.normal,
meta$cell_type.harmonized.cancer)

orig_1=as.data.frame(rna$Piece_ID)
meta=meta[rownames(orig_1),]
rna$cell_type.detailed=meta$cell_type.detailed

rna$ID=ifelse(rna$cell_type.detailed=='Tumor',paste(rna$Disease,rna$Piece_ID,rna$cell_type.detailed,sep='__'),
paste(rna$Disease, rna$cell_type.detailed,sep='__'))

var_genes <- rownames(rna)
DefaultAssay(rna)<-'SCT'

Idents(rna)<-rna$ID
aliquot.averages <- AverageExpression(rna, assays = 'SCT', slot ='data',features=var_genes)

file2write <- paste0("out/AverageExpressionSCT_SlotData_AllGenes.CellType.detailed.NormCombined.20230213.tsv")
write.table(aliquot.averages, file = file2write, quote = F, sep = "\t", row.names = T)

n_cells=as.data.frame(table(rna$ID))
write.table(n_cells, "out/Cell_count_perGroup.CellType.detailed.NormCombined.20230213.tsv",sep='\t',row.names=F,quote=F)

#########################
####ENDS HERE FOR NOW:###
#########################
