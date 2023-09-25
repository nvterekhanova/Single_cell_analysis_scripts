# Here are the scripts for generating cohort-level merged objects for pan-cancer snATAC project

## Steps:

1. ```Merge_ATAC_samples_prelim.v.1.0.R``` -- combine samples in one merged object using just small number of peaks (folder on cluster ```RDS_perCohort```).

2. ```Quantify_pancanSet_perCohort.R``` -- use pan-cancer set of peaks from this dir ```Select_peaks_for_pancan``` to quantify ATAC-coverage across samples in each cohort (folder on cluster ```Quantify_PanCanSet_perCohort```).

3. ```Process_individualRDS.100PCs.remove_Doublets.20230115.R``` -- remove doublets using tables from Scrublet and then perform normalization and dimansional reduction (folder on cluster ```Remove_doublets```).


4. ```RunChromVar.V2.R``` -- run chromVar to generate TF activity scores based on their motif accessibilites (folder on cluster ```Remove_doublets```).