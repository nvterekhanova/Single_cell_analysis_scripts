####Use counts-slot as suggested in Signac-tutorial
###Also Yige's pipeline for reference: /diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/snRNA_Processed_Data/Monocle/scripts/CellTypeVer.20201002/C3N-01200_SelfPT_SelfFib_ByCellType/
#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp

library(monocle)
library(Signac)
library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(viridis)
library(reshape)
library(reshape2)
library(plyr)
library(dplyr)


n_components=10
sample='MergedObj'
n_cells=1000

cols=c(brewer.pal(n=3, name='Dark2'),'#FF69B4','#FF0000','#00008B','#ADD8E6')
names(cols)=c('Lum progenitors','Basal progenitors','Lum mature','Her2','Basal','LumA','LumB')


annot=read.table(paste('/diskmnt/Projects/HTAN_analysis/snATAC/BR_HTAN/Signac.v.1.3/DATA/PieceID.Annotation_byCluster.Short.v4.20230503.txt',
sep=''),sep='\t',header=T)
luma=unique(annot$Piece_ID[annot$PAM50_sn %in% c('LumA')])
lumb=unique(annot$Piece_ID[annot$PAM50_sn %in% c('LumB')])
basal=unique(annot$Piece_ID[annot$PAM50_sn %in% c('Basal')])
her2=unique(annot$Piece_ID[annot$PAM50_sn %in% c('Her2')])

luma=luma[-2]
luma=c('HT088B1_S1H2','HT088B1_S1H1',luma)

all_ids=c(luma, lumb, basal, her2)


all_tab=read.table('All_PTs_monocle_objs.10comp.20230503.tsv',sep='\t',header=T)
tab_tfs=read.table('TF_score_forAll_monocleBarcodes.20230503.tsv',sep='\t',header=T)

all_tab_s=all_tab[,c('Barcode','Pseudotime','ID','Piece_ID')]

all_st_2=NULL
for (id in all_ids){

subtype='NA'
subtype=ifelse(id %in% luma, 'LumA', subtype)
subtype=ifelse(id %in% lumb, 'LumB', subtype)
subtype=ifelse(id %in% basal, 'Basal', subtype)
subtype=ifelse(id %in% her2, 'Her2', subtype)

if (subtype=='Basal'){
   cell_o='Lum progenitors'
}else{
   cell_o='Lum mature'
}

all_tab_s1=all_tab_s[all_tab_s$Piece_ID==id & all_tab_s$ID %in% c(id,cell_o),]
tab_tfs_s=tab_tfs[tab_tfs$Barcode %in% all_tab_s1$Barcode,]

res_s1=merge(tab_tfs_s,all_tab_s1,all.x=T)
res_s1$Pseudotime=as.numeric(as.character(unlist(res_s1$Pseudotime)))


all_st=NULL
for (tf in unique(res_s1$TF)){
     res_s2=res_s1[res_s1$TF==tf,]
     test=cor.test(res_s2$Pseudotime,res_s2$Score,method='pearson')
     st=cbind(tf,test$estimate,test$p.value,id)
     all_st=rbind(all_st,st)
}

all_st=as.data.frame(all_st)
colnames(all_st)=c('TF','Estimate','P_value','Piece_ID')
all_st$FDR=p.adjust(all_st$P_value, method='fdr')
all_st=all_st[order(all_st$FDR),]

all_st_2=rbind(all_st_2,all_st)
print(id)
}
write.table(all_st_2,'out/Pearson_corrs_PT_vsTFsScores.20230503.tsv',sep='\t',quote=F,row.names=F)
