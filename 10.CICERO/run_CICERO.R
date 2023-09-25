#2023-09-22: script that was used for running CICERO for ccRCC project.
system("export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp")
	       
library(BSgenome.Hsapiens.UCSC.hg38)
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(patchwork)
library(cicero)
library(EnsDb.Hsapiens.v86)

#Load RDS object
path_to_ATAC_RDS_object=''
ATAC=readRDS(path_to_ATAC_RDS_object,sep=''))

#Set the peaks assay that should be used:
DefaultAssay(ATAC)<-'peaksMACS2'

ATAC.cds <- as.cell_data_set(x = ATAC)
ATAC.cicero <- make_cicero_cds(ATAC.cds, reduced_coordinates = reducedDims(ATAC.cds)$UMAP)

#Output of the CICERO, summary statistics
#Overlap QC metrics:
#Cells per bin: 50
#Maximum shared cells bin-bin: 44
#Mean shared cells bin-bin: 0.0112756486209922
#Median shared cells bin-bin: 0

# get the chromosome sizes from the Seurat object
genome = seqlengths(BSgenome.Hsapiens.UCSC.hg38)


#convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)
chrs=paste('chr',c(1:22,'X','Y'),sep='')
genome.df=genome.df[genome.df$chr %in% chrs,]

#run cicero
conns <- run_cicero(ATAC.cicero, genomic_coords = genome.df, sample_num = 100)

write.table(conns, "out/snATAC_CICERO_conns.tsv",sep='\t',quote=FALSE,row.names=FALSE)

ccans <- generate_ccans(conns)

#From CICERO:
#"Coaccessibility cutoff used: 0.3"

write.table(ccans, "out/snATAC_CICERO_ccans.tsv",sep='\t',quote=FALSE,row.names=FALSE)

conns_1=conns[conns$coaccess>0.25,]
conns_1=conns_1[!is.na(conns_1$coaccess),]

write.table(conns_1, "out/28_ccRCC_snATAC_CICERO.0.25_cutoff.tsv",sep='\t',quote=FALSE,row.names=FALSE)

links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(ATAC) <- links

saveRDS(ATAC,snATAC_merged_obj_CICERO_added.rds,compress=T)
