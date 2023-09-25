#Select top expressed features for each cancer type, to make list of peaks for pancan-object
#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp
library(Signac)
library(Seurat)


all_peaks=NULL
for (disease in cancers){
tab=read.table(paste("../Tumor_Normal.v.20230117/features.20230117/",disease,"_features_ranking.ObjNoDoublets.TumorNormal.20230117.tsv",
sep=''),sep='\t',header=T)
all_peaks=rbind(all_peaks,tab)
print(disease)
}

st=aggregate(all_peaks$count,by=list(all_peaks$peak),FUN='sum')

#nrow(st[st$x<50,])

sum(all_peaks$count)

#for 25K:
#1862562340
#n_peaks=69158

#for 30K:
#2006388401
#n_peaks=86511


#for 35K peaks:
#2129716338
#n_peaks=103809


###Previous data freeze:
#for 35K peaks:
#1910896840
#n_peaks=103,928

#for 40K:
#2007584603
#n_peaks=121,220

#for 45K:
#2093492074
#n_peaks=138070 peaks

write.table(all_peaks, paste("Selected_features_forPanCanObject.TumorNormal.",
"Top25K.20220117.tsv",sep=''),sep='\t',row.names=F,quote=F)

###################################
##############ENDS Here for now:###
###################################
