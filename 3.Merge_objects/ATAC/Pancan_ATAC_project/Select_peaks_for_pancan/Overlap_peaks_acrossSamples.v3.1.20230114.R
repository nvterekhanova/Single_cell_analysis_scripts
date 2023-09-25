#2023-09-22: this script was used for pan-cancer snATAC paper

#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp
library(Signac)
library(Seurat)
library(future)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg38)


date='20210831'
args = commandArgs(trailingOnly=TRUE)
disease=args[1]

tab=read.table('../../ATAC_catalog.20230114.v1.txt',sep='\t',header=TRUE)

if (disease=='ccRCC'){
tab=tab[tab$Include.in.the.downstream.analysis=='TRUE' & !is.na(tab$Include.in.the.downstream.analysis=='TRUE')
 & tab$Disease.Type %in% c(disease,'PKD'),]
}else{
tab=tab[tab$Include.in.the.downstream.analysis=='TRUE' & !is.na(tab$Include.in.the.downstream.analysis=='TRUE')
 & tab$Disease.Type==disease,]
}

tab_atac=tab[tab$Data.Type=='snATAC',]
samples_atac=paste(tab_atac$Disease.Type,tab_atac$Sample.ID,sep='_')
piece_ids_atac=paste(tab_atac$Disease.Type,tab_atac$Piece_ID,sep='_')

tab_multiome=tab[tab$Data.Type=='10x_SC_Multi_ATAC_SEQ',]
samples_multiome=paste(tab_multiome$Disease.Type,tab_multiome$Sample.ID,sep='_')
piece_ids_multiome=paste(tab_multiome$Disease.Type,tab_multiome$Piece_ID,sep='_')

samples=c(samples_atac,samples_multiome)

all_peaks=NULL
for (sample in samples){
    peaks=read.table(paste("Recentered_peaks/recentered_final.filtered.v2.",sample,".tsv",
sep=""),sep='\t',header=TRUE)
    peaks$Sample=sample
    sum_score=sum(peaks$neg_log10qvalue_summit)/1000000
    peaks$norm_score=peaks$neg_log10qvalue_summit/sum_score
    all_peaks=rbind(all_peaks,peaks)
}
            
p=all_peaks

write.table(p,paste('peaks/MACS2_peaks.',disease,'_samples.BySample.',date,'.tsv',sep=''),
sep='\t',quote=FALSE,row.names=FALSE)

#p=read.table(paste('peaks/MACS2_peaks.',disease,'_samples.BySample.',date,'.tsv',sep=''),sep='\t',header=TRUE)


p1=p
rownames(p1)=c(1:nrow(p1))
recentered_p=StringToGRanges(p1$new_peak, sep = c("-", "-"))

olap=as.data.frame(findOverlaps(recentered_p,recentered_p))
olap1=olap[olap$queryHits!=olap$subjectHits,]

recentered_non_olap=p1[-olap1$queryHits,]

###Use normalized score instead (as described in Science paper), column #19
pairs=cbind(p1[olap1$queryHits,c(1,13,14,19)],olap1$queryHits,p1[olap1$subjectHits,c(1,13,14,19)],olap1$subjectHits)
colnames(pairs)=c('chr_1','st_1','en_1','score_1','row_1','chr_2','st_2','en_2','score_2','row_2')

#pairs=pairs[pairs$score_1>=pairs$score_2,]
pairs=pairs[order(-pairs$score_1),]
pairs_all=pairs

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


recentered_non_olap$row=rownames(recentered_non_olap)
recentered_final=rbind(recentered_non_olap,all_st_f)

sum_score=sum(recentered_final$norm_score)/1000000
recentered_final$norm_pancan_score=recentered_final$norm_score/sum_score

write.table(recentered_final,paste('peaks/',disease,'recentered_final.filtered.',date,'.tsv',sep=''),sep='\t',
quote=FALSE,row.names=FALSE)
write.table(recentered_non_olap,paste('peaks/',disease,'recentered_nonOverlapping.filtered.',date,'.tsv',sep=''),
sep='\t',quote=FALSE,row.names=FALSE)
write.table(all_st_f,paste('peaks/',disease,'recentered_Overlapping.filtered.',date,'.tsv',sep=''),sep='\t',
quote=FALSE,row.names=FALSE)
