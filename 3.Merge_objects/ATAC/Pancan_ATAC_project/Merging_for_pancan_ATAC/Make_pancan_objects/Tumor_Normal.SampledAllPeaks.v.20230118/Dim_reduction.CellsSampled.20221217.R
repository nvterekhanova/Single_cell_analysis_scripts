#2022-08-22: update cell type annotation for plotting.
#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp
library(Signac)
library(Seurat)
library(future)
library(Matrix)
###some parallelization-solution from the tutorial:
plan("multiprocess", workers = 40)


annot=read.table(paste('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v7.0/snATAC/Merged_objects_PanCan_peaks/',
'All_225_samples_metadata_data_freeze_v7.0.tsv',sep=''),sep='\t',header=T)
annot$cell_type.harmonized.cancer=ifelse(annot$cell_type.harmonized.cancer=='B-cells' & annot$Cancer=='MM','B-cells_MM',
annot$cell_type.harmonized.cancer)


cell_types=unique(annot$cell_type.harmonized.cancer)

###Melanocytes --- only 20 cells, so won't take them here as well ---WE WILL TAKE THEM (as was decided by the group)!!!
###Take Goblet cells, as they are normel epithelial in CRC + also add B/cells and/or Plasma (not much plasma, so just b-cells) from MM:
cell_types=unique(annot$cell_type.harmonized.cancer)
sel_cell_types=cell_types[!(cell_types %in% c('','Fibroblasts','Macrophages','T-cells','Endothelial',
'DC','Plasma','B-cells','Unknown','Pericytes','Microglia','NK','Mast','Tregs','Alveolar',
'Hepatocytes','Cholangiocytes','Neurons','Skeletal Muscle','Erythrocytes','Neurons','Skeletal Muscle','Low quality','vSMCs', 'Smooth muscle',
'Adipocytes'))]


cancers=c('OV','ccRCC','CRC','BRCA','MM','GBM','HNSCC','UCEC','PDAC','CESC','SKCM')

atac=vector(mode = "list", length = 11)
atac_old=vector(mode = "list", length = 11)

all_nnz<-NULL
n_cells=1000
for (i in 1:11){
disease=cancers[i]
atac[[i]]=readRDS(paste('../../Remove_doublets/RDS.withUMAPs.50PCs.noDoublets/',disease,
'_snATAC_Merged.PancanSet.noDoublets.20230114.rds',sep=''))


###Add_annotation
annot_sel=annot[annot$Cancer==disease,]
orig_1=as.data.frame(atac[[i]]$Piece_ID)
rownames(annot_sel)=annot_sel$V1
annot_sel=annot_sel[rownames(orig_1),]
atac[[i]]$cell_type.harmonized.cancer=annot_sel$cell_type.harmonized.cancer
atac[[i]]=subset(atac[[i]], cell_type.harmonized.cancer %in% sel_cell_types)

###Downsample object (https://satijalab.org/seurat/articles/essential_commands.html):
Idents(atac[[i]])=atac[[i]]$seurat_clusters
atac_old[[i]]=atac[[i]]
atac[[i]]=subset(x = atac[[i]], downsample = n_cells)


assay=GetAssay(atac[[i]],assay='pancan')
counts=assay@counts
#counts=counts[rownames(counts) %in% peaks,]
st=cbind(cancers[i],nnzero(counts))
all_nnz=rbind(all_nnz,st)
print(disease)
}

atac_backup=atac_old
atac=vector(mode = "list", length = 11)

all_nnz<-NULL
n_cells=1000
for (i in 1:11){
disease=cancers[i]

###Downsample object (https://satijalab.org/seurat/articles/essential_commands.html):
atac[[i]]=subset(x = atac_old[[i]], downsample = n_cells)

assay=GetAssay(atac[[i]],assay='pancan')
counts=assay@counts
st=cbind(cancers[i],nnzero(counts))
all_nnz=rbind(all_nnz,st)
print(disease)
}




all_nnz=as.data.frame(all_nnz)
sum(as.numeric(as.character(unlist(all_nnz$V2))))
#2096214430 --1K cells
#[1] 4035487701 --3K cells





combined <- merge(x = atac[[1]], y = atac[2:11])
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 'q0')

combined <- RunSVD(combined,n=200)

combined <- RunUMAP(object = combined, reduction = 'lsi', dims = 2:30,reduction.name = "umap.30PC")
combined <- RunUMAP(object = combined, reduction = 'lsi', dims = 2:50,reduction.name = "umap.50PC")
combined <- RunUMAP(object = combined, reduction = 'lsi', dims = 2:75,reduction.name = "umap.75PC")
combined <- RunUMAP(object = combined, reduction = 'lsi', dims = 2:100,reduction.name = "umap.100PC")
combined <- RunUMAP(object = combined, reduction = 'lsi', dims = 2:150,reduction.name = "umap.150PC")
combined <- RunUMAP(object = combined, reduction = 'lsi', dims = 2:200,reduction.name = "umap.200PC")

combined$Disease=gsub('(.*)_.*','\\1',combined$Piece_ID)
combined$Disease=gsub('(.*)_.*','\\1',combined$Disease)
combined$Disease=gsub('(.*)_.*','\\1',combined$Disease)
combined$Disease=ifelse(combined$Disease=='PKD','ccRCC',combined$Disease)

saveRDS(combined, "T_N_Pancan_1000cellsPerCluster_200PCs.20230118.rds.gz",compress=T)

atac=combined
DefaultAssay(atac)='pancan'

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

saveRDS(atac,"Pancan_225Samples_1000cellsPerCluster_TumorNormal.200PCs.chromVar.20230118.rds.gz",compress=T)





cols=readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v6.0/Colors_panatac_v1.0.rds')
dis_cols=cols$Cancer

####Make DimPlots:
disease='225'
p1 <- DimPlot(combined, group.by = 'Disease', pt.size = 0.1,label=T,reduction = "umap.30PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=dis_cols)
p2 <- DimPlot(combined, group.by = 'Disease', pt.size = 0.1,label=T,reduction = "umap.50PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=dis_cols)
p2_5 <- DimPlot(combined, group.by = 'Disease', pt.size = 0.1,label=T,reduction = "umap.75PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=dis_cols)
p3 <- DimPlot(combined, group.by = 'Disease', pt.size = 0.1,label=T,reduction = "umap.100PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=dis_cols)
p4 <- DimPlot(combined, group.by = 'Disease', pt.size = 0.1,label=T,reduction = "umap.150PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=dis_cols)
p5 <- DimPlot(combined, group.by = 'Disease', pt.size = 0.1,label=T,reduction = "umap.200PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=dis_cols)

pdf(paste("plots/225_Combined_snATAC.subs.diffPC_number.600CellsPerCl.20230117.pdf",sep=""),
height=7,width=8)
print(p1)
print(p2)
print(p2_5)
print(p3)
print(p4)
print(p5)
dev.off()

cell_type_cols=cols$cell_type

####DimPlots by cell type:
disease='225'
p1 <- DimPlot(combined, group.by = 'cell_type.harmonized.cancer', pt.size = 0.1,label=T,reduction = "umap.30PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=cell_type_cols)
p2 <- DimPlot(combined, group.by = 'cell_type.harmonized.cancer', pt.size = 0.1,label=T,reduction = "umap.50PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=cell_type_cols)
p2_5 <- DimPlot(combined, group.by = 'cell_type.harmonized.cancer', pt.size = 0.1,label=T,reduction = "umap.75PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=cell_type_cols)
p3 <- DimPlot(combined, group.by = 'cell_type.harmonized.cancer', pt.size = 0.1,label=T,reduction = "umap.100PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=cell_type_cols)
p4 <- DimPlot(combined, group.by = 'cell_type.harmonized.cancer', pt.size = 0.1,label=T,reduction = "umap.150PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=cell_type_cols)
p5 <- DimPlot(combined, group.by = 'cell_type.harmonized.cancer', pt.size = 0.1,label=T,reduction = "umap.200PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=cell_type_cols)

pdf(paste("plots/225_Combined_snATAC.ByCellType.subs.diffPC_number.600CellsPerCl.20230117.pdf",sep=""),
height=7,width=12)
print(p1)
print(p2)
print(p2_5)
print(p3)
print(p4)
print(p5)
dev.off()

###############################################
###2022-08-22: update cell types annotation.###
###############################################
atac=readRDS(paste('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v5.0/snATAC/Merged_objects',
'_PanCan_peaks/PanCan_merged_obj/Pancan_700cellsPerCluster_200PCs.chromVar.20220221.rds.gz',sep=''))

meta=read.table(paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/DATA/Objects.v5/',
'All_159_samples_metadata_data_freeze_v5.0.1.withDetailed.tsv',sep=''),sep='\t',header=T)

#2022-05-26: annotation needs update, because Acinar cells were re-assigned.
meta$cell_type.detailed=ifelse(!is.na(meta$cell_type.normal),meta$cell_type.normal,
meta$cell_type.harmonized.cancer)

orig_1=as.data.frame(atac$Piece_ID)
rownames(meta)=meta$X
meta=meta[rownames(orig_1),]
atac$cell_type.detailed=meta$cell_type.detailed
atac$cell_type.harmonized.cancer=meta$cell_type.harmonized.cancer

atac$Disease=ifelse(atac$Piece_ID %in% c("BRCA_HT268B1-Th1H3", "BRCA_HT029B1-S1PC", "BRCA_HT035B1-S1PA",
 "BRCA_HT1408-06","BRCA_HT141B1-S1H1", "BRCA_HT206B1-S1H4", "BRCA_HT271B1-S1H3"), "BRCA_Basal",
atac$Disease)


cols=readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v5.0/Colors_panatac_v1.0.rds')
dis_cols_ed <- c(BRCA = "#fb9a99", BRCA_Basal = "#e31a1c",CESC = "#fdbf6f",CRC = '#ff7f00',ccRCC = "#cab2d6", 
GBM = "#6a3d9a", MM = "#b15928", HNSCC = "#b2df8a",OV = "#33a02c", PDAC = "#a6cee3", UCEC = "#1f78b4" )

cell_type_cols=cols$cell_type

combined=atac
disease='159'
p1 <- DimPlot(combined, group.by = 'cell_type.harmonized.cancer', pt.size = 0.1,label=T,reduction = "umap.30PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=cell_type_cols)
p2 <- DimPlot(combined, group.by = 'cell_type.harmonized.cancer', pt.size = 0.1,label=T,reduction = "umap.50PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=cell_type_cols)
p2_5 <- DimPlot(combined, group.by = 'cell_type.harmonized.cancer', pt.size = 0.1,label=T,reduction = "umap.75PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=cell_type_cols)
p3 <- DimPlot(combined, group.by = 'cell_type.harmonized.cancer', pt.size = 0.1,label=T,reduction = "umap.100PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=cell_type_cols)
p4 <- DimPlot(combined, group.by = 'cell_type.harmonized.cancer', pt.size = 0.1,label=T,reduction = "umap.150PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=cell_type_cols)
p5 <- DimPlot(combined, group.by = 'cell_type.harmonized.cancer', pt.size = 0.1,label=T,reduction = "umap.200PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=cell_type_cols)

pdf(paste("plots/159_Combined_snATAC.ByCellType.subs.diffPC_number.700CellsPerCl.20220822.pdf",sep=""),
height=7,width=12)
print(p1)
print(p2)
print(p2_5)
print(p3)
print(p4)
print(p5)
dev.off()

###Now by Disease:
disease='159'
p1 <- DimPlot(combined, group.by = 'Disease', pt.size = 0.1,label=T,reduction = "umap.30PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=dis_cols_ed)
p2 <- DimPlot(combined, group.by = 'Disease', pt.size = 0.1,label=T,reduction = "umap.50PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=dis_cols_ed)
p2_5 <- DimPlot(combined, group.by = 'Disease', pt.size = 0.1,label=T,reduction = "umap.75PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=dis_cols_ed)
p3 <- DimPlot(combined, group.by = 'Disease', pt.size = 0.1,label=T,reduction = "umap.100PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=dis_cols_ed)
p4 <- DimPlot(combined, group.by = 'Disease', pt.size = 0.1,label=T,reduction = "umap.150PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=dis_cols_ed)
p5 <- DimPlot(combined, group.by = 'Disease', pt.size = 0.1,label=T,reduction = "umap.200PC") +
ggplot2::ggtitle(paste(disease," combined snATAC samples",sep=''))+
scale_color_manual(values=dis_cols_ed)

pdf(paste("plots/159_Combined_snATAC.subs.diffPC_number.700CellsPerCl.20220822.pdf",sep=""),
height=7,width=8)
print(p1)
print(p2)
print(p2_5)
print(p3)
print(p4)
print(p5)
dev.off()



##########################
###ENDS HERE##############
##########################
