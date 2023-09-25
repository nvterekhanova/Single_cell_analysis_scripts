#2022-01-05: Update GeneActivity normalization using same scale factor, e.g. 10000

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(RColorBrewer)

set.seed(1234)
library(plyr)
library(dplyr)
library(tibble)

library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(ggplot2)
library(future)

###some parallelization-solution from the tutorial:
plan("multiprocess", workers = 20)
options(future.globals.maxSize = 200 * 1024^3) # for 200 Gb RAM


args = commandArgs(trailingOnly=TRUE)
disease=args[1]
date='20211024'

con<-file(paste('logs/UpdateGeneActivityNormalization/',disease,'_UpdateGeneActivityNormalization.log.',date,'.txt',sep=''))
sink(con,append=TRUE)
sink(con,type='message',append=TRUE)


ATAC=readRDS(paste('out/',disease,'_snATAC_Merged.PancanSet.noDoublets.GeneActivity.',date,
'.rds',sep=''))

DefaultAssay(ATAC)='ATACGeneActivity'
options(future.globals.maxSize= 22891289600)
ATAC <- NormalizeData(
  object = ATAC,
  assay = 'ATACGeneActivity',
  normalization.method = 'LogNormalize',
  scale.factor = 10000
#  scale.factor = median(ATAC$nCount_ATACGeneActivity)
)

saveRDS(ATAC, paste('out/',disease,'_snATAC_Merged.PancanSet.noDoublets.GeneActivity.',date,
'.rds',sep=''),compress=T)
