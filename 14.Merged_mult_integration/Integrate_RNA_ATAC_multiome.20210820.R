#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(dplyr)
library(reshape)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(future)
set.seed(1234)

plan("multiprocess", workers = 40)

path_to_RNA_RDS_object=''
path_to_ATAC_RDS_object=''
samples_in_both_objs=''

#Load RNA RDS object:
rna=readRDS(path_to_RNA_RDS_object, sep=''))

#Load ATAC RDS object:
atac=readRDS(path_to_ATAC_RDS_object)

#Subset only samples from Combo (and samples that are not present in RNA RDS object)
atac_all=atac
atac=subset(atac_all, Piece_ID %in% samples_in_both_objs)

###Need to combine two objects, check if barcodes are the same in both objects:
rna_counts=GetAssay(rna, assay='RNA')
rna_counts_df=rna@assays$RNA
barcodes=colnames(rna_counts_df)
barcodes=gsub('Multiome_(.*)','\\1',barcodes)

###Since atac/rna have different barcodes, need to identify barcodes that are present in both:
orig_1=as.data.frame(atac$Piece_ID)
colnames(orig_1)[1]='Piece_ID'
orig_1$Barcode=gsub('BRCA_.*_(.*)','\\1',rownames(orig_1))
orig_1$Barcode=gsub('.*_(.*)','\\1',orig_1$Barcode)
orig_1$Multiome_barcode=paste(orig_1$Piece_ID, orig_1$Barcode,sep='_')
orig_1_subs=orig_1[orig_1$Multiome_barcode %in% barcodes,]

atac_barcodes=rownames(orig_1_subs)
atac=subset(atac, cells=atac_barcodes)

###Now change RNA-colnames to the barcode IDs used for ATAC:
orig_1_subs$Multiome_barcode_1=paste('Multiome',orig_1_subs$Multiome_barcode,sep='_')

#set order the same as in the table for ID conersion:
rna_counts_df=rna_counts_df[,orig_1_subs$Multiome_barcode_1]
colnames(rna_counts_df)=rownames(orig_1_subs)

# Now we are ready to add RNA assay. Create RNA assay and add it to the object
atac[["RNA"]] <- CreateAssayObject(
counts = rna_counts_df[,colnames(atac)])


#Need to re-run normalization for atac, since we subsetted the obj:
DefaultAssay(atac) <- "peaksMACS2"

#For this use all features, 
atac <- FindTopFeatures(atac, min.cutoff = 'q0')
atac <- RunTFIDF(atac)
atac <- RunSVD(atac)
        
atac = RunUMAP(
atac,
reduction = 'lsi',
dims = 2:30,
reduction.name = "atac.umap",
reduction.key = "atacUMAP_"
)       

DefaultAssay(atac) <- "RNA"
atac <- SCTransform(atac, vars.to.regress = c("nCount_RNA","percent.mito"),return.only.var.genes = F)
atac <- RunPCA(atac, npcs=30,verbose = FALSE)

####Try to take just embeddings:
PCA_embeddings=(Embeddings(rna, reduction = "pca"))
PCA_embeddings=PCA_embeddings[orig_1_subs$Multiome_barcode_1,]
rownames(PCA_embeddings)=rownames(orig_1_subs)
atac[["pca"]] <- CreateDimReducObject(embeddings = PCA_embeddings, assay ="SCT")


atac = RunUMAP(
atac,
reduction = 'pca',
dims = 1:30,
reduction.name = 'rna.umap',
reduction.key = 'rnaUMAP_'
)

####Now we ready to perform joint Dim-reduction:
# Joint clustering
message("[Joint] cluster cells based on joint RNA and ATAC...")
DefaultAssay(atac) = "peaksMACS2"
# build a joint neighbor graph using both assays
atac = FindMultiModalNeighbors(
object = atac,
reduction.list = list("pca", "lsi"),
dims.list = list(1:30, 2:30),
verbose = TRUE
)

# build a joint UMAP visualization
atac = RunUMAP(
object = atac,
nn.name = "weighted.nn", #weighted nearest neighbor
reduction.name = "wnn.umap",
reduction.key = "wnnUMAP_",
verbose = TRUE
)

#Also perform clustering using both modalities:
atac = FindClusters(
atac,
graph.name = "wsnn",
algorithm = 3,  # SLM
verbose = TRUE
)


#Now make plots:
p1 = DimPlot(atac, label = TRUE, repel = TRUE, reduction = "wnn.umap") +
NoLegend() +
ggtitle('Joint WNN')
p2 = DimPlot(atac, reduction = 'rna.umap', label = TRUE,
repel = TRUE, label.size = 2.5) +NoLegend() + ggtitle('RNA')
p3 = DimPlot(atac, reduction = 'atac.umap', label = TRUE,
repel = TRUE, label.size = 2.5) +NoLegend() +ggtitle('ATAC')
p = p1 + p2 + p3
pdf(paste("plots/Joint_RNA_ATAC_embeddings_Dimplots.pdf",sep=""),height=4,width=12)
print(p)
dev.off()
