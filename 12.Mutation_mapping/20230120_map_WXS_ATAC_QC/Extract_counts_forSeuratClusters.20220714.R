#2022-07-02: based on ../20220625_check_presence_WGS/Extract_counts_forTumor.20220626.R; and also using: Make_table_CNV_perGeneCluster.ATACobj.EMT.20220616.R
#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp
library(Signac)
library(Seurat)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(reshape)
library(reshape2)
library(EnsDb.Hsapiens.v86)
library(tidyverse)

#data_dir='/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/DATA/Objects.v5/'
disease='ccRCC'

###Read snATAC-annotation:
meta_final=read.table(paste('/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/DATA/Objects.v7/',
'All_225_ATAC_samples_metadata_data_freeze_v7.0.2.withDetailed.tsv',sep=''),sep='\t',header=T)

#2022-05-26: annotation needs update, because Acinar cells were re-assigned.
meta=meta_final
meta$cell_type.detailed=ifelse(!is.na(meta$cell_type.normal),meta$cell_type.normal,
meta$cell_type.harmonized.cancer)

meta$Barcode=gsub('.*\\_.*\\_(.*)','\\1',meta$X)
meta_s=meta[meta$Cancer==disease,]

all_res=NULL
piece_ids=unique(meta_s$Piece_ID)
piece_ids=piece_ids[piece_ids!='CPT704DU-M1']
#piece_ids=c('C3L-00088-T1', 'C3L-00088-T2', 'C3L-00096-T1', 'C3L-00010-T1', 'C3L-00079-T1')

if (disease=='PDAC'){
piece_ids=c('HT224P1-S1', 'HT231P1-S1H3', 'HT232P1-S1H1', 'HT242P1-S1H1', 'HT259P1-S1H1', 'HT264P1-S1H2',
 'HT270P1-S1H2', 'HT288P1-S1H4', 'HT306P1-S1H1')
}
if (disease=='CRC'){
piece_ids=c('CM618C1-T1Y2', 'CM618C2-S1Y2', 'CM663C1-S1Y1', 'CM663C1-T1Y1', 'CM1563C1-S1Y1', 'CM1563C1-T1Y1', 'CM268C1-S1', 'CM268C1-T1')
}

for (piece_id in piece_ids){
print(piece_id)

meta_s1=meta_s[meta_s$Piece_ID==piece_id,]
meta_s1=meta_s1[,c('Barcode','cell_type.detailed')]
rownames(meta_s1)=meta_s1$Barcode

muts=read.table(paste('out/',disease,'/muts_assays/res_MutAssay_',piece_id,'.txt',sep=''),sep='\t',header=T)
rownames(muts)=muts[,1]
muts=muts[,-1]
colnames(muts)=gsub('\\.','-',colnames(muts))

meta_s2=meta_s1[meta_s1$cell_type.detailed=='Tumor',]

muts_s=muts[,rownames(meta_s2)]

#Now go over each mutation and make summary:
stat1=as.data.frame(apply(muts_s,1,function(x) length(x[x==1])))
stat2=as.data.frame(apply(muts_s,1,function(x) length(x[x==2])))
colnames(stat1)='Ref'
colnames(stat2)='Alt'
stat1$mut_ID=rownames(stat1)
stat2$mut_ID=rownames(stat2)
res=merge(stat1,stat2)
res$Piece_ID=piece_id
res$Cell_type='Tumor'
res=res[order(-res$Alt),]
all_res=rbind(all_res,res)
}


all_res=as.data.frame(all_res)
all_res$Coverage=all_res$Ref+all_res$Alt
#all_res$Case=gsub('(.*)-(.*)-(.*)','\\1-\\2',all_res$Piece_ID)
write.table(all_res, paste("out/Coverage_muts_allATAC_",disease,"_samples.20230125.tsv",sep=''),sep='\t',
row.names=F,quote=F)


#Now merge with bulk:
maf=read_delim(paste('/diskmnt/Projects/snATAC_primary/06_WXS_maf/',disease,'.dnp.annotated.maf.gz',sep=''),delim='\t')

maf=as.data.frame(maf)
maf_s=maf[,c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Strand','Variant_Classification',
'Variant_Type','Reference_Allele','Tumor_Seq_Allele1','Tumor_Seq_Allele2','t_ref_count','t_alt_count',
'n_ref_count','n_alt_count','Tumor_Sample_Barcode')]

if (disease=='UCEC'){
maf_s$Piece_ID=gsub('D1_1_T','',maf_s$Tumor_Sample_Barcode)
maf_s$Piece_ID=gsub('T1Y1','T1',maf_s$Piece_ID)
maf_s$Piece_ID=gsub('S1Y1','S1',maf_s$Piece_ID)
maf_s$Piece_ID=gsub('CPT4603DU-TWIS-','',maf_s$Piece_ID)
maf_s$Case=gsub('_T','',maf_s$Piece_ID)
maf_s$Case=gsub('(.*)-H.*','\\1',maf_s$Case)
maf_s$Piece_ID=maf_s$Tumor_Sample_Barcode
}


if (disease=='CRC'){
maf_s$Piece_ID=maf_s$Tumor_Sample_Barcode
maf_s$Case=maf_s$Tumor_Sample_Barcode
maf_s$Case=gsub('(^CM[0-9]+C[0-9]+)\\-.*','\\1',maf_s$Case)
}

if (disease=='BRCA'){
   maf_s$Case=gsub('(.*)-.*','\\1',maf_s$Case)
   maf_s1=maf_s[maf_s$Case %in% c("HT035B1",atac_tab$Case[atac_tab$Disease.Type==disease]),]
}
if (disease=='UCEC'){
   maf_s$Case=gsub('(CPT.*DU)-TW.*','\\1',maf_s$Case)
}
maf_s1=maf_s[maf_s$Case %in% atac_tab$Case[atac_tab$Disease.Type==disease],]
maf_s1=maf_s1[maf_s1$Variant_Classification!='Silent',]
maf_s1$Case=gsub('CM618C.*','CM618C',maf_s1$Case)


maf_s1$mut_ID=paste(maf_s1$Hugo_Symbol, maf_s1$Chromosome, maf_s1$Start_Position,sep='_')
maf_s2=maf_s1[,c('Hugo_Symbol','Variant_Classification','t_ref_count','t_alt_count','mut_ID','Case','Piece_ID')]
colnames(maf_s2)[7]='Bulk_Piece_ID'

all_res$Case=gsub('(.*)-(.*)-(.*)','\\1-\\2',all_res$Piece_ID)
all_res_s=all_res[all_res$Coverage!=0,]

if (disease %in% c('PDAC','BRCA','UCEC','CRC')){
   all_res_s$Case=gsub('(.*)-.*','\\1',all_res_s$Case)
}
if(disease=='CRC'){
   all_res_s$Case=gsub('CM618C.*','CM618C',all_res_s$Case)
res=merge(all_res_s,maf_s2)
}

res=merge(all_res_s,maf_s2)
res=res[res$Alt>0,]



write.table(res,paste("out/Coverage_muts_allATAC_",disease,"_samples.Filtered.TumorCells.20230125.tsv",sep=''),sep='\t',quote=F,row.names=F)
############################
######ENDS MAIN CODE HERE###
############################


##########################################################################################
###Here are the scripts for extracting coverage for ATAC and Bulk (was used previously)###
##########################################################################################

bulk=read_delim('inputs/UCEC_PanCanAllNonSilent.maf',delim='\t')
bulk=as.data.frame(bulk)
bulk$Case=gsub('_T','',bulk$Tumor_Sample_Barcode)
bulk$mut_ID=paste(bulk$Hugo_Symbol, bulk$Chromosome,bulk$Start_Position,sep='_')


all_res_1=merge(all_res, bulk,all.x=T)
write.table(all_res_1, "out/Coverage_muts_allCCRCCsamples.withBulk.20220628.tsv",sep='\t',row.names=F,quote=F)

res_1_s=all_res_1
res_1_s=res_1_s[,c('Gene','t_ref_count','t_alt_count','Ref','Alt','Reference_Allele','Tumor_Seq_Allele1',
'Tumor_Seq_Allele2','mut_ID','Case','Piece_ID')]

###Calculate total coverage in ATAC Ref/Alt
stat1=aggregate(res_1_s$Ref, by=list(res_1_s$mut_ID),FUN='sum')
stat2=aggregate(res_1_s$Alt, by=list(res_1_s$mut_ID),FUN='sum')
colnames(stat1)=c('mut_ID','Ref_all_ATAC')
colnames(stat2)=c('mut_ID','Alt_all_ATAC')
stat=merge(stat1,stat2)


res_2=res_1_s[!is.na(res_1_s$Gene),]
res_3=merge(res_2,stat,all.x=T)
write.table(res_3, "out/Coverage_muts_allCCRCCsamples.withBulk.stat.20220628.tsv",sep='\t',row.names=F,quote=F)

res_1_s=res_3
res_1_s=res_1_s[res_1_s$Alt!=0,]
res_1_s=res_1_s[(res_1_s$Ref+res_1_s$Alt)>10,]
res_1_s$WGS_VAF=res_1_s$t_alt_count/(res_1_s$t_ref_count +res_1_s$t_alt_count)
res_1_s$ATAC_VAF=res_1_s$Alt/(res_1_s$Alt +res_1_s$Ref)
res_1_s$VAF_diff=res_1_s$ATAC_VAF-res_1_s$WGS_VAF
res_1_s=res_1_s[order(-res_1_s$VAF_diff),]

