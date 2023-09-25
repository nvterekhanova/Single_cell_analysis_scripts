#2023-09-22: script to generate average expression matrix per cell group. This script was used for pan-cancer snATAC project.
#2023-02-10: Fix the issue with using old harmonized.cancer
#2023-01-22: cp ../../Merged.obj.v.20221209/20221225_Peak_Coverage/20230106_test_Coverage/Peaks_coverage_test_AverageExpr.20230106.R
#2022-03-02: also extract coverage of all cell types (use cell_type.detailed).
#Peaks_coverage.v2.20220726.R
#Based on cp ../../Merged.obj.v.20220207/20220511_TumorNormalComparison/Peaks_coverage.v2.20220301.R, 
#using updated cell type annotation.
#Use slot counts:
#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp

library(BSgenome.Hsapiens.UCSC.hg38)
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
library(pheatmap)
library(viridis)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

data_dir='/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v7.0'

atac=readRDS(paste('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v7.0/snATAC/Merged_objects_PanCan_peaks/PanCan_merged_obj/Pancan_',
'225Samples_1000cellsPerCluster_TumorNormal.200PCs.chromVar.20230118.rds.gz',sep=''))

meta=read.table(paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/DATA/Objects.v7/All_225_ATAC_samples_metadata_data_freeze_v7.0.2.',
'withDetailed.tsv',sep=''),sep='\t',header=T)

meta$cell_type.detailed=ifelse(!is.na(meta$cell_type.normal),meta$cell_type.normal,
meta$cell_type.harmonized.cancer)

orig_1=as.data.frame(atac$Piece_ID)
rownames(meta)=meta$X
meta=meta[rownames(orig_1),]
atac$cell_type.detailed=meta$cell_type.detailed

###To identify top peaks -- we use separately BRCA and BRCA_Basal. For quantifying we use BRCA as one cohort, 
#as we want to combine cells from LP, BP, ML.




atac$test=ifelse(atac$cell_type.detailed=='Tumor',paste(atac$cell_type.detailed,atac$Piece_ID,sep='__'),
paste(atac$Disease, '__', atac$cell_type.detailed,sep=''))
ATAC=atac


#atac$Disease=ifelse(atac$Piece_ID %in% c("BRCA_HT268B1-Th1H3", "BRCA_HT029B1-S1PC", "BRCA_HT035B1-S1PA",
# "BRCA_HT1408-06","BRCA_HT141B1-S1H1", "BRCA_HT206B1-S1H4", "BRCA_HT271B1-S1H3",
#"BRCA_HT378B1-S1H1", "BRCA_HT378B1-S1H2", "BRCA_HT384B1-S1H1", "BRCA_HT517B1-S1H1"), "BRCA_Basal",
#atac$Disease)


cell_types=unique(ATAC$test)
Idents(ATAC)=ATAC$test

mat_access = AverageExpression(ATAC,  slot = "data")

###########################################
####How long it woud take, lets seee....###
###########################################
#quite fast!

acces=mat_access$pancan
write.table(acces, 'out/Accessibility_Piece_ID_NormCombined.DataSlot.BRCA_normal_together.20230212.tsv',sep='\t',quote=F)

orig_1=as.data.frame(table(Idents(ATAC)))
colnames(orig_1)=c('Cell_type','Cell_count')
write.table(orig_1,'out/Cell_count_perGroup.Piece_ID_NormCombined.BRCA_normal_together.20230212.tsv',sep='\t',row.names=F,quote=F)


###Now also extract data for the top 500 peaks:


daps=read.table(paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/Analysis/4.Peak_analysis/Merged.obj.v.20230114/20230121_DA_peaks/out/',
'da_peaks_oneCances_vs_Others.minPct0.1.20230126.tsv',sep=''),sep='\t',header=T)

tab=read.table('out/Accessibility_Piece_ID_NormCombined.DataSlot.20230122.tsv')



#Motifs scores for some reason doesn't correspond to the the ones we obtain by out regular script.. maybe it uses some other formulae?
###as the data are being exponentiated before calculating average -- that's why	   it doesn't correspond to our tables.
########################################
##############ENDS HERE#################
########################################
