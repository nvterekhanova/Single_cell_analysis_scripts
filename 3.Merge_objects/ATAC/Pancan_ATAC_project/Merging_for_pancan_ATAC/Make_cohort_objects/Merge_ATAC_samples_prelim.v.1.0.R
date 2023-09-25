#2023-09-22: this script was used for pan-cancer snATAC paper
#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(RColorBrewer)

require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
set.seed(1234)

library(plyr)
library(dplyr)
library(tibble)

library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(future)

###some parallelization-solution from the tutorial:
plan("multiprocess", workers = 20)
options(future.globals.maxSize = 300 * 1024^3) # for 300 Gb RAM


args = commandArgs(trailingOnly=TRUE)
disease=args[1]
date='20230114'

tab=read.table('../ATAC_catalog.20220114.v1.txt',sep='\t',header=TRUE)

if (disease=='ccRCC'){
tab=tab[tab$Include.in.the.downstream.analysis=='TRUE' & 
!is.na(tab$Include.in.the.downstream.analysis=='TRUE') & tab$Disease.Type %in% c(disease,'PKD'),]
}else{
tab=tab[tab$Include.in.the.downstream.analysis=='TRUE' & 
!is.na(tab$Include.in.the.downstream.analysis=='TRUE') & tab$Disease.Type %in% c(disease),]
}

tab_atac=tab[tab$Data.Type=='snATAC',]
samples_atac=paste(tab_atac$Disease.Type,tab_atac$Sample.ID,sep='_')
piece_ids_atac=paste(tab_atac$Disease.Type,tab_atac$Piece_ID,sep='_')

tab_multiome=tab[tab$Data.Type=='10x_SC_Multi_ATAC_SEQ',]
samples_multiome=paste(tab_multiome$Disease.Type,tab_multiome$Sample.ID,sep='_')
piece_ids_multiome=paste(tab_multiome$Disease.Type,tab_multiome$Piece_ID,sep='_')




samples=c(samples_atac,samples_multiome)
piece_ids=c(piece_ids_atac,piece_ids_multiome)

atac=vector(mode = "list", length = length(samples))
n_samples=length(samples)

if (length(samples_atac)>0){
for (i in 1:length(samples_atac)){
    atac[[i]]=readRDS(paste("/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/0.Individual_objects/snATAC/out/",disease,"/",samples[i],"/",samples[i],"_processed_atac.rds",sep=""))
    DefaultAssay(atac[[i]]) <- 'X500peaksMACS2'
    atac[[i]][['ATACGeneActivity']]<-NULL
    atac[[i]][['peaks']]<-NULL
    print (paste(i,samples[i],sep=' '))
    atac[[i]]$Piece_ID=piece_ids[i]
}
}
if (n_samples>=(length(samples_atac)+1)){
for (i in (length(samples_atac)+1):(length(samples_atac)+length(samples_multiome))){
    atac[[i]]=readRDS(paste("/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/0.Individual_objects/Multiome/out/",disease,"/",samples[i],"/",samples[i],"_processed_multiome.rds",sep=""))
    DefaultAssay(atac[[i]]) <- 'X500peaksMACS2'
    atac[[i]][['ATAC']]<-NULL
    atac[[i]][['RNA']]<-NULL
    atac[[i]][['SCT']]<-NULL
    atac[[i]][['ATACGeneActivity']]<-NULL
    print (paste(i,samples[i],sep=' '))
    atac[[i]]$Piece_ID=piece_ids[i]
}
}



#atac-assays: RNA, peaks
#X500peaksMACS2 -- make default

#from multiome: RNA, SCT,
#X500peaksMACS2 -- make default


####Make prelim merging, use same peaks:
prev_atac=readRDS('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/3.Merge_snATAC/Merge.vers.20221209/RDS_perCohort/out/UCEC_forPancan_snATAC.20221209.rds.gz')
test_peaks=StringToGRanges(rownames(prev_atac), sep = c("-", "-"))

#For testing purposes only:
peaks.use=test_peaks

matrix.counts=vector(mode = "list", length = length(samples))

for (i in 1:length(samples)){
    matrix.counts[[i]] <- FeatureMatrix(
    fragments = Fragments(atac[[i]]@assays$X500peaksMACS2),
    features = peaks.use,
    sep = c("-","-"),
    cells = colnames(atac[[i]])
    ) 
}


###add min.features=-1 (otherwise breaks)
for (i in 1:length(samples)){
atac[[i]][['peaksinters']] <- CreateChromatinAssay(counts = matrix.counts[[i]],
fragments=Fragments(atac[[i]]@assays$X500peaksMACS2),min.features=-1)
atac[[i]]$dataset=samples[i]
DefaultAssay(atac[[i]])<-'peaksinters'
###remove other assay
atac[[i]][['X500peaksMACS2']]<-NULL
}


####Merging:
combined <- merge(x = atac[[1]], y = atac[2:length(samples)], add.cell.ids = piece_ids)
saveRDS(combined, paste('out/',disease,'_forPancan_snATAC.',date,'.rds.gz',sep=''),compress=T)
