#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp
library(Signac)
library(Seurat)
library(future)
library(data.table)

cancers=c('ccRCC','BRCA','PDAC','OV','CRC','CESC','GBM','MM','UCEC','HNSCC','SKCM')

date='20210831'
date_2='20230114'

all_peaks=NULL
for (disease in cancers){
    p=read.table(paste('../',disease,'/peaks/',disease,'recentered_final.filtered.',
'20210831.tsv',sep=''),sep='\t',header=T)
    all_peaks=rbind(all_peaks,p)
}	

p1=all_peaks
rownames(p1)=c(1:nrow(p1))
recentered_p=StringToGRanges(p1$new_peak, sep = c("-", "-"))

olap=as.data.frame(findOverlaps(recentered_p,recentered_p))
olap1=olap[olap$queryHits!=olap$subjectHits,]

recentered_non_olap=p1[-olap1$queryHits,]

###Use normalized score instead (as described in Science paper), column #21, use pancan-norm score
pairs=cbind(p1[olap1$queryHits,c(1,13,14,21)],olap1$queryHits,p1[olap1$subjectHits,c(1,13,14,21)],
olap1$subjectHits)
colnames(pairs)=c('chr_1','st_1','en_1','score_1','row_1','chr_2','st_2','en_2','score_2','row_2')

###Don't remove All weaker peaks from the clusters of overlapping peaks:
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

disease='PanCan'
recentered_non_olap$row=rownames(recentered_non_olap)
recentered_final=rbind(recentered_non_olap,all_st_f)
write.table(recentered_final,paste('peaks/',disease,'recentered_final.filtered.',date_2,'.tsv',sep=''),sep='\t',
quote=FALSE,row.names=FALSE)
write.table(recentered_non_olap,paste('peaks/',disease,'recentered_nonOverlapping.filtered.',date_2,'.tsv',sep=''),
sep='\t',quote=FALSE,row.names=FALSE)
write.table(all_st_f,paste('peaks/',disease,'recentered_Overlapping.filtered.',date_2,'.tsv',sep=''),sep='\t',
quote=FALSE,row.names=FALSE)

#####ENDS HERE##############
