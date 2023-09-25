#2023-09-22: script used for merging BRCA ATAC snsamples

#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp
library(Signac)
library(Seurat)
library(future)
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(GenomeInfoDb)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
plan("multiprocess", workers = 30)
#options(future.globals.maxSize = 100 * 1024 ^ 2)

#######################################################
#############Now try the merged object, all samples:###
#######################################################

date='20211016'
n_samples=25
samples_atac=c('','','')

samples_multiome=c('','','')



####Change the score to neg_log10qvalue_summit, as it will be non-rounded
all_peaks=NULL
for (sample in samples_atac){
    peaks=read.table(paste("../../1.Create_rds/out/",sample,"/recentered_final.filtered",sample,".tsv",
sep=""),sep='\t',header=TRUE)
    peaks$Sample=sample
    sum_score=sum(peaks$neg_log10qvalue_summit)/1000000
    peaks$norm_score=peaks$neg_log10qvalue_summit/sum_score
    all_peaks=rbind(all_peaks,peaks)
}
for (sample in samples_multiome){
    peaks=read.table(paste("../../../../../Multiome/BR_HTAN/Signac.v.1.3/1.Create_rds/out/",sample,
"/recentered_final.filtered",sample,".tsv",
sep=""),sep='\t',header=TRUE)
    peaks$Sample=sample
    sum_score=sum(peaks$neg_log10qvalue_summit)/1000000
    peaks$norm_score=peaks$neg_log10qvalue_summit/sum_score
    all_peaks=rbind(all_peaks,peaks)
}
            
p=all_peaks
write.table(p,paste('peaks/MACS2_peaks.BySample.',date,'.tsv',sep=''),
sep='\t',quote=FALSE,row.names=FALSE)
p=read.table(paste('peaks/MACS2_peaks.BySample.',date,'.tsv',sep=''),sep='\t',header=TRUE)

p1=p

recentered_p=StringToGRanges(p1$new_peak, sep = c("-", "-"))

olap=as.data.frame(findOverlaps(recentered_p,recentered_p))
olap1=olap[olap$queryHits!=olap$subjectHits,]

recentered_non_olap=p1[-olap1$queryHits,]

###Use normalized score instead (as described in Science paper), column #19
pairs=cbind(p1[olap1$queryHits,c(1,13,14,19)],olap1$queryHits,p1[olap1$subjectHits,c(1,13,14,19)],
olap1$subjectHits)
colnames(pairs)=c('chr_1','st_1','en_1','score_1','row_1','chr_2','st_2','en_2','score_2','row_2')

###removing this step will allow to some more weaker peaks, as suggested here: 
#https://www.archrproject.com/bookdown/the-iterative-overlap-peak-merging-procedure.html)
#pairs=pairs[pairs$score_1>=pairs$score_2,]
pairs=pairs[order(-pairs$score_1),]
pairs_all=pairs
library(data.table)
library(doParallel)
registerDoParallel(cores=30)

p_table=p1
p_table$row=as.numeric(rownames(p_table))
p_table=as.data.table(p_table)
setkey(p_table,row)
all_st=NULL
all_st<-foreach(chr_n=c(1:22,"X","Y")) %dopar% {
chr=paste("chr",chr_n,sep='')
pairs=pairs_all[pairs_all$chr_1==chr,]
pairs=pairs[,c(4,5,9,10)]
key_list=unique(pairs$row_1)
pairs=as.data.table(pairs)
setkey(pairs,row_1)
all_st_chr=NULL
for (i in 1:length(key_list)){
    if (length(key_list)>0){
    p_add=p_table[key_list[1]]
    all_st_chr=rbind(all_st_chr,p_add)
    
    p_del=pairs[.(key_list[1])]
    del_row_1=c(p_del$row_1,p_del$row_2)
    key_list=key_list[!(key_list %in% del_row_1)] 
    }
#    print(paste(chr,length(key_list),sep=' '))
}
return(all_st_chr)
}


all_st_f=NULL
for (i in 1:24){
    all_st_1=as.data.frame(all_st[[i]])
    all_st_1=all_st_1[!duplicated(all_st_1),]
    all_st_f=rbind(all_st_f,all_st_1)
}

all_st_f=all_st_f[,1:19]
recentered_final=rbind(recentered_non_olap,all_st_f)
write.table(recentered_final,paste('peaks/recentered_final.filtered.',date,'.tsv',sep=''),sep='\t',
quote=FALSE,row.names=FALSE)
write.table(recentered_non_olap,paste('peaks/recentered_nonOverlapping.filtered.',date,'.tsv',sep=''),
sep='\t',quote=FALSE,row.names=FALSE)
write.table(all_st_f,paste('peaks/recentered_Overlapping.filtered.',date,'.tsv',sep=''),sep='\t',
quote=FALSE,row.names=FALSE)


##########################################################
####Now create FeatureMatrix with the new set of peaks:###
##########################################################
recentered_final=read.table(paste('peaks/recentered_final.filtered.',date,'.tsv',sep=''),
sep='\t',header=T)
atac=readRDS('Merged_obj_prelim.rds',sep='')

recentered_p=StringToGRanges(recentered_final$new_peak, sep = c("-", "-"))
matrix.counts <- FeatureMatrix(
    fragments = Fragments(atac@assays$peaksinters),
    features = recentered_p,
    sep = c("-","-"),
    cells = colnames(atac),
    process_n=30000
)

###next is here:
atac[['peaksMACS2']] <- CreateChromatinAssay(counts = matrix.counts,
fragments=Fragments(atac@assays$peaksinters))
DefaultAssay(atac)<-'peaksMACS2'

#atac[['X500peaksMACS2']]<-NULL
atac[['peaksinters']]<-NULL


atac <- FindTopFeatures(atac, min.cutoff = 'q0')
atac <- RunTFIDF(atac)
atac <- RunSVD(atac)

atac <- RunUMAP(object = atac, reduction = 'lsi', dims = 2:50)

p1 <- DimPlot(atac, group.by = 'Piece_ID', pt.size = 0.1) +
ggplot2::ggtitle("Combined snATAC samples")+
scale_color_manual(values=c(brewer.pal(n = 12, name = "Paired"),brewer.pal(n = 8, name = "Dark2"),
brewer.pal(n = 5, name = "Pastel1")))


pdf(paste("25_snATAC_Merged_BR.MACS2peaks_Sample.20211001.pdf",sep=""),height=6,width=8)
p1
dev.off()

#####We ran RunChromVar, and saved the object
meta=read.table(paste('/diskmnt/Projects/HTAN_analysis/snATAC/BR_HTAN/Signac.v.1.2/Annotation/19_BR_snATAC.',
'Merged.CellTypeAnnotation.v1.20210719.tsv',sep=''),sep='\t',header=T)
ATAC=atac
orig_1=as.data.frame(ATAC$Piece_ID)
colnames(orig_1)='Piece_ID'
orig_1$Barcode=gsub('.*_(.*)','\\1',row.names(orig_1))
orig_1$Barcode=paste(orig_1$Piece_ID,orig_1$Barcode,sep='__')
meta$Barcode=gsub('.*_(.*)','\\1',row.names(meta))
meta$Barcode=paste(meta$Piece_ID,meta$Barcode,sep='__')
rownames(meta)=meta$Barcode
meta=meta[orig_1$Barcode,]
ATAC$cell_type=meta$cell_type


p1 <- DimPlot(ATAC, group.by = 'cell_type', pt.size = 0.1) +
ggplot2::ggtitle("25 BRCA snATAC samples")+
scale_color_manual(values=c(brewer.pal(n = 10, name = "Paired")))

pdf(paste("plots/25_snATAC_Merged_BR.MACS2peaks_byCellTypePrelim.20211016.pdf",sep=""),height=6,width=8)
p1
dev.off()

ATAC <- FindNeighbors(object = ATAC, reduction = 'lsi', dims = 2:50)

###this helps:
options(future.globals.maxSize= 891289600)
ATAC <- FindClusters(object = ATAC, verbose = FALSE, algorithm = 3)

p1 <- DimPlot(ATAC, group.by = 'seurat_clusters', pt.size = 0.1, label=T) +
ggplot2::ggtitle("25 BRCA snATAC samples")

pdf(paste("plots/25_snATAC_Merged_BR.MACS2peaks_bySeuratClusters.20211016.pdf",sep=""),height=6,width=8)
p1
dev.off()

###Do Manual cell_type annotation
ATAC$cell_type_manual='Tumor'
ATAC$cell_type_manual=ifelse(ATAC$seurat_clusters %in% c(9,19),'mCAF',ATAC$cell_type_manual)
ATAC$cell_type_manual=ifelse(ATAC$seurat_clusters %in% c(35),'vCAF',ATAC$cell_type_manual)
ATAC$cell_type_manual=ifelse(ATAC$seurat_clusters %in% c(31),'DC',ATAC$cell_type_manual)
ATAC$cell_type_manual=ifelse(ATAC$seurat_clusters %in% c(0),'Macrophage',ATAC$cell_type_manual)
ATAC$cell_type_manual=ifelse(ATAC$seurat_clusters %in% c(1,10),'T_NK',ATAC$cell_type_manual)
ATAC$cell_type_manual=ifelse(ATAC$seurat_clusters %in% c(26),'B',ATAC$cell_type_manual)
ATAC$cell_type_manual=ifelse(ATAC$seurat_clusters %in% c(12),'Plasma',ATAC$cell_type_manual)
ATAC$cell_type_manual=ifelse(ATAC$seurat_clusters %in% c(25),'Endothelial',ATAC$cell_type_manual)
ATAC$cell_type_manual=ifelse(ATAC$seurat_clusters %in% c(37),'Mast',ATAC$cell_type_manual)
#ATAC$cell_type_manual=ifelse(ATAC$seurat_clusters %in% c(),'',ATAC$cell_type_manual)

p1 <- DimPlot(ATAC, group.by = 'cell_type_manual', pt.size = 0.1, label=T) +
ggplot2::ggtitle("25 BRCA snATAC samples")+
scale_color_manual(values=c(brewer.pal(n = 10, name = "Paired")))

pdf(paste("plots/Merged_snATAC.MACS2peaks_byCellTypeManual.pdf",sep=""),height=6,width=8)
p1
dev.off()

saveRDS(ATAC, "Merged_snATAC_obj.rds.gz",compress=T)

####Add GeneActivity:
ATAC=readRDS("Merged_snATAC_obj.rds.gz")
#Caclulte GeneActivitty for the merged ATAAC object

#genome(ATAC)<-'hg38'
DefaultAssay(ATAC)<-'peaksMACS2'
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
genome(annotations) <- "NA"

# change to UCSC style since the data was mapped to hg19; should be this order:
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(ATAC) <- annotations

gene.activities <- GeneActivity(ATAC)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
ATAC[['ATACGeneActivity']] <- CreateAssayObject(counts = gene.activities)

options(future.globals.maxSize= 22891289600)
ATAC <- NormalizeData(
  object = ATAC,
  assay = 'ATACGeneActivity',
  normalization.method = 'LogNormalize',
  scale.factor = median(ATAC$nCount_ATACGeneActivity)
)

saveRDS(ATAC, "Merged_snATAC_obj.GeneActivity_added.rds.gz",compress=T)