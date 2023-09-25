#2023-09-22: this script was used for pan-cancer snATAC paper
#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp

library(Signac)
library(Seurat)
library(GenomeInfoDb)

require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
set.seed(1234)

library(plyr)
library(dplyr)
library(tibble)
library(future)

###some parallelization-solution from the tutorial:
plan("multiprocess", workers = 20)
options(future.globals.maxSize = 200 * 1024^3) # for 200 Gb RAM


args = commandArgs(trailingOnly=TRUE)
disease=args[1]
date='20230114'

peaks=read.table('../Select_peaks.20230114/PanCan/peaks/PanCanrecentered_final.filtered.20230114.tsv',
sep='\t',header=T)
atac=readRDS(paste('../RDS_perCohort/out/',
disease,'_forPancan_snATAC.',date,'.rds.gz',sep=''))

peaks.use=StringToGRanges(peaks$new_peak, sep = c("-", "-"))

    matrix.counts <- FeatureMatrix(
    fragments = Fragments(atac@assays$peaksinters),
    features = peaks.use,
    sep = c("-","-"),
    cells = colnames(atac),
    process_n=30000
    ) 


###add min.features=-1 (otherwise breaks)
atac[['pancan']] <- CreateChromatinAssay(counts = matrix.counts,
fragments=Fragments(atac@assays$peaksinters),min.features=-1)
DefaultAssay(atac)<-'pancan'

#Extract features presence across cells into separate files
features=FindTopFeatures(atac[['pancan']], min.cutoff = 'q0')
features_s=as.data.frame(features@meta.features)
features_s=features_s[order(-features_s$count),]
features_s$peak=rownames(features_s)
features_s$Disease=disease

write.table(features_s, paste("features/",disease,"_features_ranking.tsv",sep=''),sep='\t',quote=F,row.names=F)

saveRDS(atac, paste("RDS/",disease,"_snATAC_Merged.PancanSet.",date,".rds",sep=''))

#####ENDS HERE