system("export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp")

###For some samples (with many cells >6K) python-package used in chromVar doesn't work properly; Need to use this 2 commands:
###export OMP_NUM_THREADS=1
###export USE_SIMPLE_THREADED_LEVEL3=1
###export OPENBLAS_NUM_THREADS=1

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



ATAC=readRDS(paste('/diskmnt/Projects/HTAN_analysis/snATAC/BR_HTAN/Signac.v.1.3/3.Merge_snATAC/',
'Merge.SelectPeaks.20211016/25_BR_snATAC.selectedPeaks.chromvar.MotifsAdded.CellTyped.v1.20211017.rds.gz',
sep=''))

####Use cell_type_manual for cell_type Annotation

ATAC$cell_type_manual=ifelse(ATAC$seurat_clusters==13,"Lum mature",ATAC$cell_type_manual)
ATAC$cell_type_manual=ifelse(ATAC$seurat_clusters==15,"Basal progenitors",ATAC$cell_type_manual)
ATAC$cell_type_manual=ifelse(ATAC$seurat_clusters==30,"Lum progenitors",ATAC$cell_type_manual)
atac=ATAC
ATAC=subset(atac, cell_type_manual %in% c('Tumor', "Lum mature", "Basal progenitors", "Lum progenitors"))


ATAC$test=paste(ATAC$cell_type_manual, ATAC$Piece_ID,sep='__')
ATAC$test=ifelse(ATAC$cell_type_manual %in% c("Lum mature", "Basal progenitors", "Lum progenitors"), ATAC$cell_type_manual, ATAC$test)


DefaultAssay(ATAC) <- 'chromvar'

chromv= GetAssayData(object = ATAC)
jaspar=read.table('/diskmnt/Projects/HTAN_analysis/snATAC/Signac/CCRCC/JASPAR2020_motifs.txt',sep='\t',
header=TRUE)

mtx0=chromv
res=merge(mtx0,jaspar,by=0,all.x=TRUE)
rownames(res)=res$motif.name
res=res[,-1]
res=res[,1:(ncol(res)-2)]


ann_col0 = data.frame('cell_types' = ATAC$test)
ann_col0$cel_2=ann_col0$cell_types
ann_col0=ann_col0[order(ann_col0$cell_types),]
ann_col1=data.frame("cell_types"=ann_col0$cell_types)
rownames(ann_col1)=rownames(ann_col0)

res=res[,rownames(ann_col0)]

final_wilcoxon_stat=NULL
final_wilcoxon_stat_all=NULL

cell_types=unique(ATAC$test)
for (cell_type in cell_types){
    if(length(rownames(ann_col1)[ann_col1$cell_types==cell_type])>=50){
    res_1=res[,colnames(res) %in% rownames(ann_col1)[ann_col1$cell_types==cell_type]]
    all_wilcoxon_stat=NULL
        for (motif in 1:nrow(res)){
           mean_score=mean(as.numeric(as.character(unlist(res_1[motif,]))),na.rm=TRUE)
           stat=cbind(cell_type,rownames(res)[motif],mean_score)
           all_wilcoxon_stat=rbind(all_wilcoxon_stat,stat)
        }
        all_wilcoxon_stat=as.data.frame(all_wilcoxon_stat)
        all_wilcoxon_stat$mean_score=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score)))
        all_wilcoxon_stat=all_wilcoxon_stat[order(-all_wilcoxon_stat$mean_score),]
        final_wilcoxon_stat=rbind(final_wilcoxon_stat,all_wilcoxon_stat)
        print(cell_type)
}
}

final_wilcoxon_stat=as.data.frame(final_wilcoxon_stat)
colnames(final_wilcoxon_stat)[2]='TF_Name'

write.table(final_wilcoxon_stat,paste("out/Motif_score_perCellGroup.for6a.20230509.tsv"),
quote=FALSE,sep="\t",row.names=FALSE)



####################ENDS HERE FOR NOW##############