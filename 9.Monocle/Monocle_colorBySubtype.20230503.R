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

all_tab=NULL
for (id in all_ids){
print(id)

subtype='NA'
subtype=ifelse(id %in% luma, 'LumA', subtype)
subtype=ifelse(id %in% lumb, 'LumB', subtype)
subtype=ifelse(id %in% basal, 'Basal', subtype)
subtype=ifelse(id %in% her2, 'Her2', subtype)

atac=readRDS(paste('../Monocle_RDS/1000_MergedObj_',id,'_',subtype,'_10_monocle.rds',sep=''))


p1=plot_cell_trajectory(atac, color_by = "Pseudotime", cell_size = 1) +
scale_color_viridis(option="B")

pdf(paste('plots.pt/',id,'_',subtype,'_monocle.byPT.20230503.pdf',sep=''),width=6,height=5,useDingbats=F)
print(p1)
dev.off()

p2=plot_cell_trajectory(atac, color_by = "ID", cell_size = 1) +
scale_color_manual(values=cols)

pdf(paste('plots.subtype/',id,'_',subtype,'_monocle.bySubtype.20230503.pdf',sep=''),width=6,height=5,useDingbats=F)
print(p2)
dev.off()

tab1=cbind(sampleNames(atac),atac$Pseudotime,atac$State,atac$test)
tab1=as.data.frame(tab1)
colnames(tab1)=c('Barcode','Pseudotime','State','ID')
tab1$Piece_ID=id
all_tab=rbind(all_tab,tab1)
}
write.table(all_tab,'All_PTs_monocle_objs.10comp.20230503.tsv',sep='\t',quote=F,row.names=F)


###Extract TF scores:
atac=readRDS(paste('/diskmnt/Projects/HTAN_analysis/snATAC/BR_HTAN/Signac.v.1.3/3.Merge_snATAC/Merge.SelectPeaks.20211016/25_BR_snATAC.',
'selectedPeaks.chromvar.MotifsAdded.CellTyped.GeneActivity.v2.20211018.rds.gz',sep=''))
DefaultAssay(atac)='chromvar'
chromv=GetAssayData(object=atac)

jaspar=read.table('/diskmnt/Projects/HTAN_analysis/snATAC/Signac/CCRCC/JASPAR2020_motifs.txt',sep='\t',
header=TRUE)
mtx0=chromv
res=merge(mtx0,jaspar,by=0,all.x=TRUE)
rownames(res)=res$motif.name
res=res[,-1]
res=res[,1:(ncol(res)-2)]


res_1=res[,colnames(res) %in% unique(all_tab$Barcode)]
res_2=melt(as.matrix(res_1))
colnames(res_2)=c('TF','Barcode','Score')
write.table(res_2,'TF_score_forAll_monocleBarcodes.20230503.tsv',sep='\t',quote=F,row.names=F)

#length(unique(all_tab$Barcode)) #N=25,478 unique barcodes
