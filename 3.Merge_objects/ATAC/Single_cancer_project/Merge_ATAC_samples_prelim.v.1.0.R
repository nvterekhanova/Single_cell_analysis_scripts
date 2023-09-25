#2023-09-22: script used for merging BRCA ATAC snsamples
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(RColorBrewer)

set.seed(1234)

library(plyr)
library(dplyr)
library(tibble)
library(reshape)
library(plyr)

library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(future)

###some parallelization-solution from the tutorial:
plan("multiprocess", workers = 30)
options(future.globals.maxSize = 300 * 1024^3) # for 300 Gb RAM
date='20211016'

n_samples=25
samples_atac=c('','','')

samples_multiome=c('','','')


tab=read.table('Case_id_table.20211016.txt',sep='\t',header=T)
rownames(tab)=tab$Sample
tab1=tab[samples_atac,]
piece_ids_atac=tab1$Case

tab1=tab[samples_multiome,]
piece_ids_multiome=tab1$Case

samples=c(samples_atac,samples_multiome)
piece_ids=c(piece_ids_atac,piece_ids_multiome)

atac=vector(mode = "list", length = length(samples))

for (i in 1:length(samples_atac)){
    atac[[i]]=readRDS(paste("../../2.Remove_doublets/out/",samples[i],
"_processed_atac.noDoublets.rds",sep=""))
    DefaultAssay(atac[[i]]) <- 'X500peaksMACS2'
    atac[[i]][['ATACGeneActivity']]<-NULL
    atac[[i]][['peaks']]<-NULL
    print (paste(i,samples[i],sep=' '))
    atac[[i]]$Piece_ID=piece_ids[i]
}
for (i in (length(samples_atac)+1):(length(samples_atac)+length(samples_multiome))){
    atac[[i]]=readRDS(paste("../../../../../Multiome/BR_HTAN/Signac.v.1.3/2.Remove_doublets/out/",
samples[i],"_processed_multiome.noDoublets.rds",sep=""))
    DefaultAssay(atac[[i]]) <- 'X500peaksMACS2'
    atac[[i]][['ATACGeneActivity']]<-NULL
    atac[[i]][['ATAC']]<-NULL
    atac[[i]][['RNA']]<-NULL
    atac[[i]][['SCT']]<-NULL
    print (paste(i,samples[i],sep=' '))
    atac[[i]]$Piece_ID=piece_ids[i]
}


combined.peaks <- UnifyPeaks(object.list = atac, mode = "reduce")
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks
peaks.use=sample(combined.peaks, size = 5000, replace = FALSE)


#We don't filter cells like in the tutorial, because we use already filtered matrices. And all cells would pass those filters in the tutorial.

matrix.counts=vector(mode = "list", length = length(samples))


for (i in 1:length(samples)){
    matrix.counts[[i]] <- FeatureMatrix(
    fragments = Fragments(atac[[i]]@assays$X500peaksMACS2),
    features = peaks.use,
    sep = c("-","-"),
    cells = colnames(atac[[i]]),
    process_n=30000
    ) 
}


for (i in 1:length(samples)){
atac[[i]][['peaksinters']] <- CreateChromatinAssay(counts = matrix.counts[[i]],
fragments=Fragments(atac[[i]]@assays$X500peaksMACS2))
atac[[i]]$dataset=samples[i]
DefaultAssay(atac[[i]])<-'peaksinters'
###remove other assay
#atac[[i]][['X500peaksMACS2']]<-NULL
}


####Merging:
combined <- merge(x = atac[[1]], y = atac[2:length(samples)], add.cell.ids = samples)
saveRDS(combined, paste(n_samples,'_mergedObj_snATAC.',date,'.rds',sep=''))
