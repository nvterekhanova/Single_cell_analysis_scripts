####Important!!! --to limit number of cores:
library(BiocParallel)
register(SerialParam())
register(MulticoreParam(20)) #number of cores to use - otherwise it crushes

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

cancers=c('OV','ccRCC','CRC','BRCA','MM','GBM','HNSCC','UCEC','PDAC','CESC','SKCM')

for (can in cancers){

print(can)
atac=readRDS(paste('RDS.withUMAPs.50PCs.noDoublets/',can,'_snATAC_Merged.PancanSet.noDoublets.20230114.rds',sep=''))

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


saveRDS(atac, paste("RDS.withUMAPs.50PCs.noDoublets.chromVar/",can,"_snATAC_Merged.PancanSet.noDoublets.chromVar.20230118.rds",sep=''),
compress=T)
}