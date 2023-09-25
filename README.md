# Single_cell_analysis_scripts

## Scripts used for the analysis of single cell data: snATAC, snRNA, and multiome.

   * Scripts folder contains general scripts used for data processing and analysis.

## General pipeline scripts for snATAC data reviewed on 2023-09-22

### Software versions:

  * Cellranger-atac v.2

  * Cellranger-arc v.2

  * Cellranger v.6.0.2

  * Signac v.1.3.0 (should be at least 1.2.1 to be able to work with cellranger-arc v.2 results)

  * Seurat v.4.0.3

##Consider using the ones from the latest versions.

## Steps to process snATAC-data (10X genomics)

### I. Run Cellranger-atac:

   * Organise symlimlinks for .fq files in the separate folder "FASTQ"

   * Download cellranger-atac and reference files: https://support.10xgenomics.com/single-cell-atac/software/downloads/latest

   * With this params it will occupy ~0.25 of one of katmai machines for less than a day (exact time depends on sample)

   * Check here for the help in interpretation of cellranger-atac web_summary.html report: https://assets.ctfassets.net/an68im79xiti/Cts31zFxXFXVwJ1lzU3Pc/fe66343ffd3039de73ecee6a1a6f5b7b/CG000202_TechnicalNote_InterpretingCellRangerATACWebSummaryFiles_Rev-A.pdf


##### Command to run Cellranger-atac (v.2.0.0), example for BR sample; reference is the same as for cellranger-arc v.2:

```/diskmnt/Software/cellranger-atac-2.0.0/cellranger-atac count --id 1408-06 --fastqs /diskmnt/Projects/HTAN_analysis_2/Cellranger-atac/BR_HTAN/Cellranger_atac_v.2.0/preprocessing/1408-06/FASTQ/ --reference /diskmnt/Datasets/Reference/CellRanger-ATAC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 --localcores=12 --localmem=100```



### II.a. Create RDS object using atac_pipeline_v3.2.R based on the Signac package (v.1.3.0 was used for testing)

  * Install latest Signac using instructions here: https://satijalab.org/signac/articles/install.html.

   + I recommend to install it in the separaate conda env. Based on the feedback it can be complicated - ask others for the help if needed (ss_atac_seq/single_cell channels).

  * Check the Signac tutorial based on which the pipeline was created: https://satijalab.org/signac/articles/pbmc_vignette.html. Also check for the changes: /diskmnt/Projects/HTAN_analysis/snATAC/BR_HTAN/Signac.v.1.3/Scripts/snATAC/README

  * Install via conda MACS2 and provide path "/home/nvterekhanova/anaconda3/envs/r-environment/bin/macs2"

  * Copy "hg38.chrom.sizes.txt" to the working dir.

##### Example run of the pipeline (the params below usually work fine for the majrity of samples; also they were used for pancan-atac study using samples across 10 cancer types):

```
sample="sample_id"

cellranger_output_directory=/diskmnt/Projects/HTAN_analysis_2/Cellranger-atac/BR_HTAN/Cellranger_atac_v.2.0

macs2_path=/home/nvterekhanova/anaconda3/envs/r-environment/bin/macs2

Rscript atac_pipeline_v3.3.R -s $sample --atac_data $cellranger_output_directory --macs2_path $macs2_path --prf_min 1000 --pct_min 15 --ns_max 5 --atac_pc_first 2

--sample #sample_name

--atac_data #path to data folder (e.g. cellranger output's raw matrices folder)

--prf_min #peak_region_fragments minimum value for filtering of cells

--prf_max #peak_region_fragments maximum value for filtering of cells

--pct_min #pct_reads_in_peaks_minimum value for filtering, percentage of reads in peaks minimum

--bl_ratio #blacklist_ratio_minimum value for filtering

--ns_max #nucleosome_signal_maximum value for filtering

--tss #tss_enrichment_minimum value for filtering

--pc_num #number of principal components to use

--atac_pc_first #first principal components to use (should be 1 or 2)
```



### II.b. Create RDS object using atac_pipeline_v1.2.R based on the Signac package (v.1.3.0 was used for testing)

  * Check the Signac tutorial based on which the pipeline was created: https://satijalab.org/signac/articles/pbmc_multiomic.html; also check for the changes: /diskmnt/Projects/HTAN_analysis/snATAC/BR_HTAN/Signac.v.1.3/Scripts/Multiome/README

##### Example run of the pipeline (the params below usually work fine for the majrity of samples; also they were used for pancan-atac study using samples across 10 cancer types)

```
sample="sample_id"

cellranger_output_directory=/diskmnt/Projects/HTAN_analysis_2/Cellranger-arc/rerun_v2.0

macs2_path=/home/nvterekhanova/anaconda3/envs/r-environment/bin/macs2

Rscript multiome_pipeline_v1.1.R -s $sample --atac_data $cellranger_output_directory --macs2_path $macs2_path --prf_min 1000 --pct_min 15 --ns_max 5 --atac_pc_first 2
```

  * Use parameters as for the atac_pipeline + parameters used for rna_pipeline (Seurat_v3_analysis_auto_v0.6.R)



### III. Annotate RDS object using already annotated sn/scRNA RDS-object

  * It will transfer the labels from the snRNA.rds meta.data field "cell_type" (change accordingly) to the snATAC.rds object, and will save the annotated ATAC.rds and the plots with annotation into the out/ dir.

##### Example run of the script:

```Rscript Integrate_withRNA.R $atac_sample.rds $rna_sample.rds```

  * This step is not needed for merging, and can be done later.



### IV. Remove doublets using scrublet (https://github.com/swolock/scrublet)

  * This step can be done before or after merging; but in case if after -- you will need to re-do normalization and Dim. reduction steps.



### V. Merge snATAC objects:

  * Create merged object from individual snATAC-objects. This will also find overlapping set of peaks across samples and will re-calculate coverage using this new set of peaks for each cell.

  * Script to use (run in the ineractive mode): ```Merge_ATAC_samples_auto.v.3.0.R```

   + Edit input samples before running it



### VI. Run ChromVar to calculate motif scores for all cells:

  * Script to use, this new version will also save motif positions mapped to peaks to the RDS object -- save the object (run script in the interactive mode): ```RunChromVar.v2.R```.
