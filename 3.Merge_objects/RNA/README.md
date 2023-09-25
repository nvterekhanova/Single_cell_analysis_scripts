# Merge snRNA cohort objects, and make PanCan object.

  * Also can be used for merging of individual snRNA samples.

## Analysis steps.

1. ```Input_merged_RNA_RDS``` -- Make symlinks here to cohort level snRNA objects that you need to merge.

2. ```Merge_snRNA_samples.20230118.R``` -- Merging of cohort level object, and make PanCan object.

3. ```Top_markers_expression.20211215.R``` -- Calculate average gene expression per cell group of top variable genes.
