####Important!!! --to limit number of cores:
library(BiocParallel)
register(SerialParam())
register(MulticoreParam(30)) #number of cores to use - otherwise it crushes on some clusters

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)

library(reshape)
library(plyr)
library(dplyr)

library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

library(ggplot2)
library(RColorBrewer)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
set.seed(1234)

#Set the path to ATAC RDS file:
path_to_ATAC_RDS_file=''

atac=readRDS(path_to_ATAC_RDS_file)

DefaultAssay(atac)='peaksMACS2'
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

saveRDS(atac,"Merged_snATAC.selectedPeaks.chromvar.MotifsAdded.rds")