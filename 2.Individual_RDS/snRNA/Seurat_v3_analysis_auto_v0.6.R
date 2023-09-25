#!/usr/bin/env Rscript --vanilla

# load required libraries
library(optparse)
library(Seurat)
library(dplyr)
library(Matrix)
library(RColorBrewer)
library(ggplot2)

# create user options
option_list = list(
  make_option(c("-i", "--input"),
              type="character",
              default=NULL,
              help="path to data folder (e.g. cellranger output's raw matrices folder)",
              metavar="character"),
  make_option(c("--pre_filter"),
              type="integer",
              default=300,
              help="min number of reads per cell to prefilter",
              metavar="integer"),
  make_option(c("--nfeature_min"),
              type="integer",
              default=200,
              help="nFeature_RNA min value for filtering",
              metavar="integer"),
  make_option(c("--nfeature_max"),
              type="integer",
              default=10000,
              help="nFeature_RNA max value for filtering",
              metavar="integer"),
  make_option(c("--ncount_min"),
              type="integer",
              default=1000,
              help="nCount_RNA min value for filtering",
              metavar="integer"),
  make_option(c("--ncount_max"),
              type="integer",
              default=80000,
              help="nCount_RNA max value for filtering",
              metavar="integer"),
  make_option(c("--mito_max"),
              type="double",
              default=0.1,
              help="maximum allowed mitochondrial fraction",
              metavar="double"),
  make_option(c("-o", "--output"),
              type="character",
              default="./",
              help="output folder path",
              metavar="character"),
  make_option(c("-s","--sample_id"),
              type="character",
              default="single_cell_study",
              help="Name of your sample",
              metavar="character"),
  make_option(c("--pc_num"),
              type="integer",
              default=30,
              help="number of principal components to use",
              metavar="integer")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# complain if there's no data
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Path to data is required (--input).n", call.=FALSE)
}

# read in initial arguments
sample_id <- opt$sample_id
out_path <- opt$output
matrix_dir = opt$input

# make output dir if it doesn't exist
dir.create(out_path)

# get direct paths to data
barcode.path <- paste0(matrix_dir, "/barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "/features.tsv.gz")
matrix.path <- paste0(matrix_dir, "/matrix.mtx.gz")

# read in matrix
input <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE,stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE,stringsAsFactors = FALSE)
colnames(input) = barcode.names$V1
rownames(input) = feature.names$V2

# pre-filter and create initial seurat object
bc_sums <- colSums(input)
bg_id <- names(bc_sums[bc_sums >= opt$pre_filter])
panc = CreateSeuratObject(counts = input[,bg_id],project=opt$sample_id,min.cells = 0)

### QC
# get percent mitochondrial content
mito.genes <- grep(pattern = "^MT-", x = feature.names$V2, value=F)
percent.mito <- Matrix::colSums((GetAssayData(panc,slot="counts"))[mito.genes, ])/Matrix::colSums(GetAssayData(panc,slot="counts"))

# plot pre-filter metadata
panc$percent.mito<-percent.mito
pdf(paste(out_path,"/QC_in_sample_",sample_id, ".pdf", sep=""), width=15, height=9)
VlnPlot(object = panc, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

# plot metadata associations
pdf(paste0(out_path,"/FeatureScatter_in_sample_",sample_id,".pdf",sep=""),width=12,height=7)
plot1 <- FeatureScatter(object = panc, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(object = panc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

# filter step
panc<-subset(x = panc, subset = nFeature_RNA > opt$nfeature_min & nFeature_RNA < opt$nfeature_max & nCount_RNA > opt$ncount_min & nCount_RNA < opt$ncount_max & percent.mito<opt$mito_max)
panc %>% dim

# plot post-filter metadata
pdf(paste(out_path,"/After_QC_in_sample_",sample_id, ".pdf", sep=""), width=15, height=9)
VlnPlot(object = panc, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), ncol = 3)
dev.off()

# Run the standard workflow for visualization and clustering
panc <- SCTransform(panc, vars.to.regress = c("nCount_RNA","percent.mito"),return.only.var.genes = F)
panc <- RunPCA(panc, npcs = opt$pc_num, verbose = FALSE)

# t-SNE and Clustering
panc <- RunUMAP(panc, reduction = "pca", dims = 1:opt$pc_num)
panc <- FindNeighbors(panc, reduction = "pca", dims = 1:opt$pc_num)
panc <- FindClusters(panc, resolution = 0.5)

# plot the clusters
pdf(paste0(out_path,"/DimPlot_",sample_id,".pdf"),useDingbats=FALSE)
DimPlot(object = panc, reduction = "umap",label=TRUE,label.size=6)
dev.off()

# save object so far
saveRDS(panc,file = paste(out_path,"/",sample_id, "_processed.rds", sep=""))
