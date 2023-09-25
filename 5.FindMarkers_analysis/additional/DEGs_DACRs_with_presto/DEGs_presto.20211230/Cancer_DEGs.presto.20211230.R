#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp
library(Signac)
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(ggplot2)
RhpcBLASctl::blas_set_num_threads(50)
library(future)
plan("multiprocess", workers = 5)
options(future.globals.maxSize = 50 * 1024^3)

path_to_RNA_RDS=''
rna=readRDS(path_to_RNA_RDS,sep=''))
rna_all=rna
rna=subset(rna_all, cell_type.harmonized.cancer=='Tumor')

rna$Piece_ID=paste(rna$Disease,rna$Piece_ID, sep='_')
basal_samples_ids=''
rna$Disease=ifelse(rna$Piece_ID %in% basal_samples_ids, "BRCA_Basal",
rna$Disease)
dis_cols_ed=c("BRCA"= "#E9967A","BRCA_Basal"="#C70039", "CESC"="#FFFF00", "CRC"="#FF8C00", "GBM"="#6A3D9A",
"HNSCC"="#FF69B4", "MM"="#A65628", "OV"="#57C785", "PDAC"="#80B1D3","UCEC"="#1F78B4", "ccRCC"="#C0C0C0")


###Now try running presto! using normalized counts:
DefaultAssay(rna)<-'SCT'
markers_disease <- presto:::wilcoxauc.Seurat(X = rna, group_by = 'Disease', assay = 'data',
seurat_assay = 'SCT')
markers_disease_sel=markers_disease[markers_disease$padj<0.05 & markers_disease$logFC>0.1,]

write.table(markers_disease, 'out/DEGs_byDisease_All.presto.20211230.tsv',sep='\t',quote=F, row.names=F)
write.table(markers_disease_sel, 'out/DEGs_byDisease_padj.0.05.logFC.0.1.presto.20211230.tsv',sep='\t',
quote=F, row.names=F)


####Alternative way of calculation DEGs with FindMarkers by Seurat:

cancers=unique(rna$Disease)
Idents(rna)=rna$Disease
DefaultAssay(rna)='SCT'
all_degs=NULL

for (disease in cancers){
degs <- FindMarkers(
  object = rna,
  ident.1 = disease,
#  ident.2='',
  only.pos = FALSE,
  min.pct = 0.1,
  min.diff.pct=0,
  assay='SCT',
  logfc.threshold=0
)

degs$Disease=disease
degs$Gene=rownames(degs)
all_degs=rbind(all_degs,degs)
print(disease)
}
write.table(all_degs, paste("out/degs_oneCances_vs_Others.minPct0.1.part1.20211229.tsv",sep=""),
sep="\t",quote=FALSE,row.names=FALSE)
write.table(all_degs, paste("out/degs_oneCances_vs_Others.minPct0.1.20211229.tsv",sep=""),
sep="\t",quote=FALSE,row.names=FALSE)
