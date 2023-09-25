#2023-09-22: script used for subtype DAMs calculation for BRCA project.
system("export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp")

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


###Cluters:
#13 - Lum mature
#15 - Basal prog
#30 - Lum prog
ATAC=readRDS(paste('/diskmnt/Projects/HTAN_analysis/snATAC/BR_HTAN/Signac.v.1.3/3.Merge_snATAC/',
'Merge.SelectPeaks.20211016/25_BR_snATAC.selectedPeaks.chromvar.MotifsAdded.CellTyped.v1.20211017.rds.gz',
sep=''))

####Use cell_type_manual for cell_type Annotation

ATAC$cell_type_manual=ifelse(ATAC$seurat_clusters==13,"Lum mature",ATAC$cell_type_manual)
ATAC$cell_type_manual=ifelse(ATAC$seurat_clusters==15,"Basal progenitors",ATAC$cell_type_manual)
ATAC$cell_type_manual=ifelse(ATAC$seurat_clusters==30,"Lum progenitors",ATAC$cell_type_manual)
atac=ATAC
ATAC=subset(atac, cell_type_manual %in% c('Tumor') & Piece_ID!='HT384B1_M1')
#ATAC=subset(atac, cell_type_manual %in% c('Tumor','Lum mature','Basal progenitors','Lum progenitors'))


annot=read.table('../PieceID.Annotation_byCluster.Short.v3.20211202.txt',sep='\t',header=T)
lum=unique(annot$Piece_ID[annot$PAM50_sn %in% c('LumA','LumB')])
lum=gsub('HT088B1','HT088B1_S1H1',lum)
lum=c(lum,'HT088B1_S1H2')
basal=unique(annot$Piece_ID[annot$PAM50_sn %in% c('Basal')])
her2=unique(annot$Piece_ID[annot$PAM50_sn %in% c('Her2')])

ATAC$test=ATAC$Piece_ID
ATAC$test=ifelse(ATAC$Piece_ID %in% lum, 'Lum',ATAC$test)
ATAC$test=ifelse(ATAC$Piece_ID %in% basal, 'Basal',ATAC$test)
ATAC$test=ifelse(ATAC$Piece_ID %in% her2, 'Her2',ATAC$test)

#lum=paste(lum, 'Tumor',sep='_')
#basal=paste(basal, 'Tumor',sep='_')



####Run the comparisons:
#cols=c("Basal"="#FF0000",'Her2'='#FF69B4','LumA'='#00008B','LumB'='#ADD8E6','Normal'='#008000','NA'='#D3D3D3')

DefaultAssay(ATAC) <- 'chromvar'

chromv= GetAssayData(object = ATAC)
#cell_types=unique(as.numeric(as.character(unlist(ATAC$test))))
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

###First Basal Tumors vs Lum progenitors
cell_types=her2
for (cell_t1 in cell_types){
    print (cell_t1)

    cell_t2='Lum'
    res_1=res[,colnames(res) %in% rownames(ann_col1)[ann_col1$cell_types==cell_t1]]
    res_2=res[,colnames(res) %in% rownames(ann_col1)[ann_col1$cell_types==cell_t2]]
    all_wilcoxon_stat=NULL
        for (motif in 1:nrow(res)){
           mean_score1=mean(as.numeric(as.character(unlist(res_1[motif,]))),na.rm=TRUE)
           mean_score2=mean(as.numeric(as.character(unlist(res_2[motif,]))),na.rm=TRUE)
           w_test=wilcox.test(as.numeric(as.character(unlist(res_1[motif,]))),
as.numeric(as.character(unlist(res_2[motif,]))))
           stat=cbind(cell_t1,rownames(res)[motif],mean_score1,mean_score2,w_test$p.value)
           all_wilcoxon_stat=rbind(all_wilcoxon_stat,stat)
        }
        all_wilcoxon_stat=as.data.frame(all_wilcoxon_stat)
        all_wilcoxon_stat$V5=as.numeric(as.character(unlist(all_wilcoxon_stat$V5)))
        all_wilcoxon_stat=all_wilcoxon_stat[order(all_wilcoxon_stat$V5),]
        all_wilcoxon_stat$mean_score1=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score1)))
        all_wilcoxon_stat$mean_score2=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score2)))
        all_wilcoxon_stat$diff=all_wilcoxon_stat$mean_score1-all_wilcoxon_stat$mean_score2
        all_wilcoxon_stat=all_wilcoxon_stat[order(-all_wilcoxon_stat$diff),]

colnames(all_wilcoxon_stat)[1:2]=c('cell_t1','TF')
colnames(all_wilcoxon_stat)[5]=c('pvalue')
all_wilcoxon_stat$FDR=p.adjust(all_wilcoxon_stat$pvalue, method="fdr")


final_wilcoxon_stat=rbind(final_wilcoxon_stat,all_wilcoxon_stat)
}

final_wilcoxon_stat=as.data.frame(final_wilcoxon_stat)
colnames(final_wilcoxon_stat)[1:2]=c('cell_t1','TF_Name')
final_wilcoxon_stat$cell_t2=cell_t2

write.table(final_wilcoxon_stat,paste("out/Score_difference.Her2_vsLum.20211203.tsv",
sep=""),quote=FALSE,sep="\t",row.names=FALSE)




######################################
### Lum/Basa/Her2 vs All others ######
######################################
cell_types=c('Basal','Her2','Lum')
for (cell_t1 in cell_types){
    print (cell_t1)

    res_1=res[,colnames(res) %in% rownames(ann_col1)[ann_col1$cell_types==cell_t1]]
    res_2=res[,colnames(res) %in% rownames(ann_col1)[ann_col1$cell_types!=cell_t1]]
    all_wilcoxon_stat=NULL
        for (motif in 1:nrow(res)){
           mean_score1=mean(as.numeric(as.character(unlist(res_1[motif,]))),na.rm=TRUE)
           mean_score2=mean(as.numeric(as.character(unlist(res_2[motif,]))),na.rm=TRUE)
           w_test=wilcox.test(as.numeric(as.character(unlist(res_1[motif,]))),
as.numeric(as.character(unlist(res_2[motif,]))))
           stat=cbind(cell_t1,rownames(res)[motif],mean_score1,mean_score2,w_test$p.value)
           all_wilcoxon_stat=rbind(all_wilcoxon_stat,stat)
        }
        all_wilcoxon_stat=as.data.frame(all_wilcoxon_stat)
        all_wilcoxon_stat$V5=as.numeric(as.character(unlist(all_wilcoxon_stat$V5)))
        all_wilcoxon_stat=all_wilcoxon_stat[order(all_wilcoxon_stat$V5),]
        all_wilcoxon_stat$mean_score1=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score1)))
        all_wilcoxon_stat$mean_score2=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score2)))
        all_wilcoxon_stat$diff=all_wilcoxon_stat$mean_score1-all_wilcoxon_stat$mean_score2
        all_wilcoxon_stat=all_wilcoxon_stat[order(-all_wilcoxon_stat$diff),]

colnames(all_wilcoxon_stat)[1:2]=c('cell_t1','TF')
colnames(all_wilcoxon_stat)[5]=c('pvalue')
all_wilcoxon_stat$FDR=p.adjust(all_wilcoxon_stat$pvalue, method="fdr")


final_wilcoxon_stat=rbind(final_wilcoxon_stat,all_wilcoxon_stat)
}

final_wilcoxon_stat=as.data.frame(final_wilcoxon_stat)
colnames(final_wilcoxon_stat)[1:2]=c('cell_t1','TF_Name')
final_wilcoxon_stat$cell_t2='Other'

write.table(final_wilcoxon_stat,paste("out/Score_difference.EachSubt_vs_Others.20220105.tsv",
sep=""),quote=FALSE,sep="\t",row.names=FALSE)
