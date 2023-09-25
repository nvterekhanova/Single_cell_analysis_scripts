#2023-09-22: script used for merging RNA cohort-level objects to generate pan-cancer RNA object
###https://satijalab.org/seurat/articles/sctransform_vignette.html
#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp
library(Signac)
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(ggplot2)
RhpcBLASctl::blas_set_num_threads(60)

date='20230118'

cancers=c("SKCM","HNSCC","CESC","OV","UCEC","ccRCC","BRCA","CRC","MM","PDAC","GBM")
###Also read cell type annotation:

path_to_RNA_meta_data=''
annot=read.table(path_to_RNA_meta_data, sep='\t', header=T)

rna=vector(mode = "list", length = length(cancers))

for (i in 1:length(cancers)){
     disease=cancers[i]
    rna[[i]]=readRDS(paste("../Input_merged_RNA_RDS/",disease,"_merged_obj.rds",sep=""))
    DefaultAssay(rna[[i]]) <- 'RNA'
    rna[[i]][['SCT']]<-NULL
    rna[[i]]$Disease=cancers[i]

    print(disease)
}

annot$cell_type.harmonized.cancer=ifelse(annot$cell_type.harmonized.cancer=='B-cells' & annot$Cancer=='MM','B-cells_MM',
annot$cell_type.harmonized.cancer)

cell_types=unique(annot$cell_type.harmonized.cancer)
sel_cell_types=cell_types[!(cell_types %in% c('','Fibroblasts','Macrophages','T-cells','Endothelial',
'DC','Plasma','B-cells','Unknown','Pericytes','Microglia','NK','Mast','Tregs','Alveolar',
'Hepatocytes','Cholangiocytes','Neurons','Skeletal Muscle','Erythrocytes','Neurons','Skeletal Muscle','Low quality','vSMCs',
'Adipocytes','Immune','Neutrophils','Monocytes','Pre-B-cells','Smooth muscle'))]


rna_backup=rna

for (i in 1:length(cancers)){
    disease=cancers[i]
    ###Add_annotation
    annot_sel=annot[annot$Cancer==disease,]
    orig_1=as.data.frame(rna[[i]]$Piece_ID)
    rownames(annot_sel)=annot_sel[,1]
    annot_sel=annot_sel[rownames(orig_1),]
    rna[[i]]$cell_type.harmonized.cancer=annot_sel$cell_type.harmonized.cancer    
#    rna[[i]]=subset(rna[[i]],cell_type.harmonized.cancer %in% sel_cell_types)
    print(disease)
}

###Check barcode matching: 
for (i in 1:11){
 print(nrow(annot[annot$Cancer==cancers[i],]))
 }
for (i in 1:11){
 print(sum(table(rna[[i]]$cell_type.harmonized.cancer)))
 }
#they match 100%!

for (i in 1:length(cancers)){
    disease=cancers[i]
    rna[[i]]=subset(rna[[i]],cell_type.harmonized.cancer %in% sel_cell_types)
    print(disease)
}



###in snATAC we used "add.cell.ids =samples" --apply same later to snATAC as well
obj = merge(x=rna[[1]],y=rna[2:length(cancers)])

DefaultAssay(obj) = "RNA"

#Seems that methd chice is not critical, so just use the one used for individual cohorts:
obj <- SCTransform(obj,vars.to.regress = c("nCount_RNA","percent.mito"), 
return.only.var.genes = F,conserve.memory=T)

###Try different method:https://satijalab.org/seurat/articles/sctransform_vignette.html
#obj <- SCTransform(obj, method = "glmGamPoi", vars.to.regress = c("nCount_RNA","percent.mito"), 
#return.only.var.genes = F,conserve.memory=T)

### try different number of PCs:
obj <- RunPCA(obj, npcs = 200)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30,reduction.name = "umap.30PC")
obj <- RunUMAP(obj, reduction = "pca", dims = 1:50,reduction.name = "umap.50PC")
obj <- RunUMAP(obj, reduction = "pca", dims = 1:75,reduction.name = "umap.75PC")
obj <- RunUMAP(obj, reduction = "pca", dims = 1:100,reduction.name = "umap.100PC")
obj <- RunUMAP(obj, reduction = "pca", dims = 1:150,reduction.name = "umap.150PC")
obj <- RunUMAP(obj, reduction = "pca", dims = 1:200,reduction.name = "umap.200PC")

obj <- FindNeighbors(obj, reduction = "pca", dims = 1:50, force.recalc = T)
obj <- FindClusters(obj, resolution = 0.5)

saveRDS(obj, "Merged_206_snRNA.SCTrandformed.3KVarFeatures.TumorNormal.v7.20230118.rds")


cols=readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v6.0/Colors_panatac_v1.0.rds')

obj$Piece_ID=paste(obj$Disease,obj$Piece_ID, sep='_')
obj$Disease=ifelse(obj$Piece_ID %in% c("BRCA_HT268B1-Th1H3", "BRCA_HT029B1-S1PC", "BRCA_HT035B1-S1PA",
 "BRCA_HT1408-06","BRCA_HT141B1-S1H1", "BRCA_HT206B1-S1H4", "BRCA_HT271B1-S1H3"), "BRCA_Basal",
obj$Disease)
dis_cols_ed=cols$Cancer

####Make DimPlots:
disease='206'
p1 <- DimPlot(obj, group.by = 'Disease', pt.size = 0.1,label=T,reduction = "umap.30PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=dis_cols_ed)
p2 <- DimPlot(obj, group.by = 'Disease', pt.size = 0.1,label=T,reduction = "umap.50PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=dis_cols_ed)
p2_5 <- DimPlot(obj, group.by = 'Disease', pt.size = 0.1,label=T,reduction = "umap.75PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=dis_cols_ed)
p3 <- DimPlot(obj, group.by = 'Disease', pt.size = 0.1,label=T,reduction = "umap.100PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=dis_cols_ed)
p4 <- DimPlot(obj, group.by = 'Disease', pt.size = 0.1,label=T,reduction = "umap.150PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=dis_cols_ed)
p5 <- DimPlot(obj, group.by = 'Disease', pt.size = 0.1,label=T,reduction = "umap.200PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=dis_cols_ed)

pdf(paste("plots/206_Combined_snATAC.subs.diffPC_number.pdf",sep=""),height=7,width=8)
print(p1)
print(p2)
print(p2_5)
print(p3)
print(p4)
print(p5)
dev.off()

######now color by cell type:
disease='206'
p1 <- DimPlot(obj, group.by = 'cell_type.harmonized.cancer', pt.size = 0.1,label=T,reduction = "umap.30PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=cols$cell_type)
p2 <- DimPlot(obj, group.by = 'cell_type.harmonized.cancer', pt.size = 0.1,label=T,reduction = "umap.50PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=cols$cell_type)
p2_5 <- DimPlot(obj, group.by = 'cell_type.harmonized.cancer', pt.size = 0.1,label=T,reduction = "umap.75PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=cols$cell_type)
p3 <- DimPlot(obj, group.by = 'cell_type.harmonized.cancer', pt.size = 0.1,label=T,reduction = "umap.100PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=cols$cell_type)
p4 <- DimPlot(obj, group.by = 'cell_type.harmonized.cancer', pt.size = 0.1,label=T,reduction = "umap.150PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=cols$cell_type)
p5 <- DimPlot(obj, group.by = 'cell_type.harmonized.cancer', pt.size = 0.1,label=T,reduction = "umap.200PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=cols$cell_type)

pdf(paste("plots/206_Combined_snATAC.byCellType.diffPC_number.pdf",sep=""),height=7,width=13)
print(p1)
print(p2)
print(p2_5)
print(p3)
print(p4)
print(p5)
dev.off()


###to check cell specific markers across different cell types:
###plot markers:
Plasma = c("SDC1","IGHG1","IGHG3","IGHG4")
B_cell = c("CD19","MS4A1","CD79A","CD79B","CD83","CD86")
Monocyte = c("LYZ","CD14","S100A8")
Myeloid = c("FCGR3A", "MS4A7", "IFITM3","CSF1R","ITGAM","CD163","CSF3R","S100A12","PTPRC","EPCAM","CD68","ANPEP",
"ITGAX","CD14","CD33","CD4")
DC = c("FCER1A","CD40","PTPRC","IL3RA","LILRA4","CD1C","CLEC7A","CLEC4C","ITGAE","ITGAM","ITGAX","LY75","CD8A")
CD8_T = c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B","EVL","CD69","CCL5","ARHGDIB","CD7","IL2RB","RAC2","PTPRC")
CD4_T = c("CD3D", "CD3E", "CD3G", "CD6","IL7R", "LDHB", "NOSIP", "CD4","TBX21","PTPRC")
T_naive = c("CD4","CD8A","CD8B","CCR7","SELL","LEF1")
T_reg = c("CD4","IL2RA","FOXP3","CTLA4","TIGIT","TNFRSF4","LAG3","PDCD1","HAVCR2")
NK = c("GNLY", "TYROBP", "HOPX", "FCGR3A","CXCL2","NCR1","PTPRC","NCAM1","KLRC1","PFN1","GZMA","GZMB","GZMH",
"IFNG","PRF1")
Erythrocyte = c("HBD", "GYPA", "HBA1", "HBA2", "CA1","HBB","BRSK1")
Endothelial = c("VWF", "PECAM1", "FLT4", "FLT1", "FLT3", "KDR","PLVAP","ANGPT2","TRIM24","ACTA2")
Epithelial = c("EPCAM","AMBP","MUC1")
Mast = c("TPSB2","TPSD1","TPSAB1","ENPP3","KIT")
Stroma = c("EPCAM","VIM","CD44","KRT8","C4BPA","COL1A1","COL1A2")
Fibroblast = c("TIMP1","FN1","POSTN","ACTA2","BST2","LY6D","COL6A1","SLC20A1","COL6A2","KRT16","CD9","S100A4",
"EMP1","LRRC8A","EPCAM","PDPN","ITGB1","PDGFRA","THY1")
Tuft = c("AVIL","BMX","AOC1","KCTD12","SH2D6","LIPG","DKK3","PLG","CORO2B")
Epithelial = c("KRT19","KRT8","KRT18","KRT17","KRT7","KRT5","KRT6A","KRT14","TACSTD2","ANXA2","S100A10",
"S100A11","S100A16","TPM1","TFF1","S100A6","AGR2","C19orf33","PDPN","ACTA2","PDGFRA","GSTP1")

ATAC=obj1
x1=DotPlot(ATAC, features=Myeloid)+ggtitle(paste('Myeloid',sep=''))+ RotatedAxis()
x2=DotPlot(ATAC, features=Endothelial)+ggtitle(paste('Endothelial',sep=''))+ RotatedAxis()
x3=DotPlot(ATAC, features=CD8_T)+ggtitle(paste('CD8_T',sep=''))+ RotatedAxis()
x4=DotPlot(ATAC, features=CD4_T)+ggtitle(paste('CD4_T',sep=''))+ RotatedAxis()
x5=DotPlot(ATAC, features=Fibroblast)+ggtitle(paste('Fibroblast',sep=''))+ RotatedAxis()
x6=DotPlot(ATAC, features=T_reg)+ggtitle(paste('T_reg',sep=''))+ RotatedAxis()
x7=DotPlot(ATAC, features=Stroma)+ggtitle(paste('Stroma',sep=''))+ RotatedAxis()
x8=DotPlot(ATAC, features=Epithelial)+ggtitle(paste('Epithelial',sep=''))+ RotatedAxis()
x9=DotPlot(ATAC, features=B_cell)+ggtitle(paste('B_cell',sep=''))+ RotatedAxis()
x10=DotPlot(ATAC, features=Plasma)+ggtitle(paste('Plasma',sep=''))+ RotatedAxis()
x11=DotPlot(ATAC, features=NK)+ggtitle(paste('NK',sep=''))+ RotatedAxis()
#x9=DotPlot(ATAC, features=)+ggtitle(paste('',sep=''))+ RotatedAxis()
#x9=DotPlot(ATAC, features=)+ggtitle(paste('',sep=''))+ RotatedAxis()

pdf(paste('plots/Markers_expr.pdf',sep=''),width=10,height=7)
print(x1)
print(x2)
print(x3)
print(x4)
print(x5)
print(x6)
print(x7)
print(x8)
print(x9)
print(x10)
print(x11)
dev.off()