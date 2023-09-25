#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp
library(Signac)
library(Seurat)
library(sctransform)
library(glmGamPoi)
RhpcBLASctl::blas_set_num_threads(50)

var_genes <- VariableFeatures(obj)
DefaultAssay(obj)<-'SCT'

obj$ID=paste(obj$cell_type.harmonized.cancer, obj$Disease, sep='_')
Idents(obj)<-obj$ID
aliquot.averages <- AverageExpression(obj, assays = 'SCT', slot ='data',features=var_genes)

file2write <- paste0("out/AverageExpressionSCT_SlotData_3KariableGenes.20211215.tsv")
write.table(aliquot.averages, file = file2write, quote = F, sep = "\t", row.names = T)
