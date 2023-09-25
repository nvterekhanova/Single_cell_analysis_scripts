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

path_to_ATAC_RDS=''
atac=readRDS(path_to_ATAC_RDS,sep=''))

ATAC=atac
atac=subset(ATAC, cell_type.harmonized.cancer=='Tumor')
####First test presto performance:
Idents(atac)=atac$Disease


#####Try using normalized counts-matrix:
assay=GetAssay(atac,assay='pancan_s')
counts=assay@counts
peak_counts=colSums(counts)
counts=t(counts)
counts_norm=counts/peak_counts
counts_norm=counts_norm*10000
counts_norm=t(counts_norm)

###Replace_counts slot (for testing purposes), probably better to create separate assay:
#atac <- SetAssayData(atac, slot = "counts", new.data = counts_norm)


atac[['pancan_snorm']] <- CreateChromatinAssay(counts = counts_norm,
fragments=Fragments(atac@assays$pancan_s),min.features=-1)

###Now try running presto! using normalized counts:
DefaultAssay(atac)<-'pancan_snorm'
markers_disease <- presto:::wilcoxauc.Seurat(X = atac, group_by = 'Disease', assay = 'counts', 
seurat_assay = 'pancan_snorm')
markers_disease_sel=markers_disease[markers_disease$padj<0.05 & markers_disease$logFC>0.1,]

write.table(markers_disease,"out/markers_disease.20211030.tsv",sep='\t',quote=F,row.names=F)
write.table(markers_disease_sel,"out/markers_disease.padj0.05.logFC.0.1.20211030.tsv",sep='\t',
quote=F,row.names=F)

####Annotate those peaks:
peaks_1=StringToGRanges(markers_disease_sel$feature, sep = c("-", "-"))
 ###Now annotate peaks:
peakAnno <- annotatePeak(peaks_1, tssRegion=c(-1000, 100),
                         TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db")
                         
anno=as.data.frame(peakAnno)
peaks=markers_disease_sel
peaks$Gene=anno$SYMBOL
peaks$Type=anno$annotation
peaks$geneId=anno$geneId
peaks$peak_distanceToTSS=anno$distanceToTSS
write.table(peaks,"out/markers_disease.padj0.05.logFC.0.1.Annotated.20211030.tsv",sep='\t',
quote=F,row.names=F)