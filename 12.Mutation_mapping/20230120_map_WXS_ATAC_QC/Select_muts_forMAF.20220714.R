#2022-07-14: we need to check just some samples (6 with low CNV-fractions), but maybe check them all.
#2022-07-01: based on Overlap_withPromoters.20220625.R.
#2022-06-25: evaluate how many inteeresting instancess we have.
#use for annotation code from: /diskmnt/Projects/Users/nvterekhanova/Germline_analysis/3.Annotate_SNPS
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

disease='UCEC'

data_dir='/diskmnt/Projects/Users/nvterekhanova/PanCan_snATAC/DATA/Objects.v5/'

atac_tab=read.table(paste('ATAC_catalog.BAM_paths.20230120.v1.txt',sep=''),sep='\t',header=T)


maf=read_delim(paste('/diskmnt/Projects/snATAC_primary/06_WXS_maf/',disease,'.dnp.annotated.maf.gz',sep=''),delim='\t')

maf=as.data.frame(maf)
maf_s=maf[,c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Strand','Variant_Classification',
'Variant_Type','Reference_Allele','Tumor_Seq_Allele1','Tumor_Seq_Allele2','t_ref_count','t_alt_count',
'n_ref_count','n_alt_count','Tumor_Sample_Barcode')]

maf_s$Piece_ID=gsub('D1_1_T','',maf_s$Tumor_Sample_Barcode)
maf_s$Piece_ID=gsub('T1Y1','T1',maf_s$Piece_ID)
maf_s$Piece_ID=gsub('S1Y1','S1',maf_s$Piece_ID)
maf_s$Piece_ID=gsub('CPT4603DU-TWIS-','',maf_s$Piece_ID)
maf_s$Case=gsub('_T','',maf_s$Piece_ID)
maf_s$Case=gsub('(.*)-H.*','\\1',maf_s$Case)
if (disease=='CRC'){
maf_s$Case=gsub('(.*)-.*-.*','\\1',maf_s$Case)
}

if (disease=='BRCA'){
   maf_s$Case=gsub('(.*)-.*','\\1',maf_s$Case)
   maf_s1=maf_s[maf_s$Case %in% c("HT035B1",atac_tab$Case[atac_tab$Disease.Type==disease]),]
}
if (disease=='UCEC'){
   maf_s$Case=gsub('(CPT.*DU)-TW.*','\\1',maf_s$Case)
}
if (disease=='CRC'){
   maf_s$Case=gsub('(.*)-.*','\\1',maf_s$Case)
}

maf_s1=maf_s[maf_s$Case %in% atac_tab$Case[atac_tab$Disease.Type==disease],]
maf_s1=maf_s1[maf_s1$Variant_Classification!='Silent',]


maf_test=maf[maf$Tumor_Sample_Barcode %in% maf_s1$Tumor_Sample_Barcode & 
maf$Start_Position %in% maf_s1$Start_Position,]

write.table(maf_test, paste('inputs/',disease,'/',disease,'_PanCanAllNonSilent.maf',sep=''),sep='\t',
row.names=F,quote=F)



atac_tab_s=atac_tab[atac_tab$Disease.Type==disease & atac_tab$Sample.Type %in% c('Tumor','Met','Normal'),]

for (piece_id in c(atac_tab_s$Piece_ID)){

path=atac_tab_s$Data.path[atac_tab_s$Piece_ID==piece_id]

path_id=paste('inputs/',disease,'/BAMs/',piece_id,'.possorted_bam.bam',sep='')
system(paste("ln -s ", path, " ",  path_id,sep=''))

path_bai=paste(path,'.bai',sep='')
path_bai_id=paste('inputs/',disease,'/BAMs/',piece_id,'.possorted_bam.bam.bai',sep='')
system(paste("ln -s ", path_bai, " ",  path_bai_id,sep=''))
print(piece_id)
}


for (piece_id in atac_tab_s$Piece_ID){
dir.create(paste('out/',disease,'/muts.mapped/',piece_id,sep=''))
}