# Pipeline for processing snRNA data by Dan.

  * Creates indidual RDS object after QCs filtering of cells.

## Notes:

### Usage: Seurat_v3_analysis_auto_v0.5.R [options]


### Options:

```
-i CHARACTER, --input=CHARACTER path to data folder (e.g. cellranger output's raw matrices folder)

-pf INTEGER, --pre_filter=INTEGER
	min number of reads per cell to prefilter

-fmin INTEGER, --nfeature_min=INTEGER
	nFeature_RNA min value for filtering

-fmax INTEGER, --nfeature_max=INTEGER
	nFeature_RNA max value for filtering

-cmin INTEGER, --ncount_min=INTEGER
	nCount_RNA min value for filtering

-cmax INTEGER, --ncount_max=INTEGER
	nCount_RNA max value for filtering

-m DOUBLE, --mito_max=DOUBLE
	maximum allowed mitochondrial fraction

-o CHARACTER, --output=CHARACTER
	output folder path

-s CHARACTER, --sample_id=CHARACTER
	Name of your sample

-sc INTEGER, --scale_factor=INTEGER
	scale factor for normalization

-v INTEGER, --variable_feature=INTEGER
	number of variable features to consider

-pc INTEGER, --pc_num=INTEGER
	number of principal components to use

-h, --help
	Show this help message and exit
```