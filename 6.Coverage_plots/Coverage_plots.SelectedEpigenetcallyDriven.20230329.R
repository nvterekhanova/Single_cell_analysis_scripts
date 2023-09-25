#2023-09-22: this script was used for pan-cancer snATAC paper.
#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp
library(Signac)
library(Seurat)
library(presto)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(reshape)
library(reshape2)
library(EnsDb.Hsapiens.v86)
library(tidyverse)



path_to_ATAC_RDS_object=''
path_to_ATAC_meta_data=''

path_to_RNA_RDS_object=''
path_to_RNA_meta_data=''
atac=readRDS(path_to_ATAC_RDS_object,sep=''))


meta=read.table(path_to_ATAC_meta_data, sep='\t', header=T)

meta$cell_type.detailed=ifelse(!is.na(meta$cell_type.normal),meta$cell_type.normal,
meta$cell_type.harmonized.cancer)

orig_1=as.data.frame(atac$Piece_ID)
rownames(meta)=meta$X
meta=meta[rownames(orig_1),]
atac$cell_type.detailed=meta$cell_type.detailed


basal_samples_ids='
atac$Disease=ifelse(atac$Piece_ID %in% basal_samples_ids, "BRCA_Basal", atac$Disease)


####Annotate those peaks:
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
genome(annotations) <- "NA"

# change to UCSC style
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(atac) <- annotations


sel_cell_o=c('Normal squamous cells','Ductal-like2','Distal Stem Cells','B-cells','Luminal mature','Luminal progenitor',
'Secretory Endometrial epithelial cells','OPC','Proximal Tubule','Melanocytes')

#B-cells

ATAC=atac
atac=subset(ATAC, cell_type.detailed %in% c('Tumor'))

Idents(atac)=atac$Disease
Idents(atac)=factor(Idents(atac),levels=c('HNSCC','CESC','PDAC','CRC','MM','BRCA','BRCA_Basal','OV',
'UCEC','GBM','ccRCC','SKCM'))

###2022-06-03: update annotation:
data_dir='/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v7.0'

rna=readRDS(paste(data_dir, path_to_RNA_RDS_object)

meta=read.table(path_to_RNA_meta_data,sep='\t',header=T)
rownames(meta)=meta$X
meta$cell_type.detailed=ifelse(!is.na(meta$cell_type.normal),meta$cell_type.normal,
meta$cell_type.harmonized.cancer)

orig_1=as.data.frame(rna$Piece_ID)
meta=meta[rownames(orig_1),]
rna$cell_type.detailed=meta$cell_type.detailed

rna$Disease=ifelse(rna$Piece_ID %in% basal_samples_ids, "BRCA_Basal", rna$Disease)


RNA=rna
rna=subset(RNA, cell_type.detailed %in% c('Tumor'))

Idents(rna)=rna$Disease
Idents(rna)=factor(Idents(rna),levels=c('HNSCC','CESC','PDAC','CRC','MM','BRCA','BRCA_Basal','OV',
'UCEC','GBM','ccRCC','SKCM'))

#2022-08-02: run across all selected peaks:

annot=Annotation(atac)
sq_genes=c('KRT6A','KRT5','TNS4','KRT15','KRT17','S100A9','ARHGEF4','SYCP2','PLD1','FABP5')
ov_ucec_genes=c('FGFR1','MAGI1','TPRG1','SLC40A1','WDR45B','TMTC1','PAX8','LYPD6B','MEIS1','WFDC2')
genes_selected=c(sq_genes,ov_ucec_genes)

for (gene in genes_selected){
###extract gene coords
region=as.data.frame(LookupGeneCoords(atac, gene, assay = NULL))
strand=as.character(unique(strand(annot[annot$gene_name==gene,])))[1]
chr=region$seqnames
if (strand=='+'){
st=region$start
}else if (strand=='-'){
st=region$end
}
en=region$end
new_st=as.numeric(st)-2000
new_en=as.numeric(st)+2000
new_peak=paste(chr,new_st,new_en,sep='-')
print(paste(gene,' ',strand,sep=''))
expr_plot <- ExpressionPlot(
  object = rna,
  features = gene,
  assay = "SCT",
  slot="data"
)
cov_plot=CoveragePlot(
  object = atac,
  region = new_peak,
  annotation = FALSE,
  peaks = FALSE,
  links=FALSE
)
gene_plot <- AnnotationPlot(
  object = atac,
  region = new_peak
)
p=CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  widths=c(3.5,1)
)

pdf(paste("plots/CoverageByDisease_",peak,"_",gene,".pdf",sep=""),width=4.5,height=5.5,useDingbats=FALSE)
print(p)
dev.off()
}
