# Analysis of pseudotime trajectories using Monocle2.

  * For each Monocle2 run we will have a separate folder, here ```Monocle.run.v.20211211```.

## Notes:

1. ```Monocle_DDR.MergedObj.v5.20211211.R``` -- Script for running Monocle2 analysis.

  * Use counts-slot as suggested in Signac tutorial

  * Here we used 10 components for param ```max_components``` everywhere. Usually trying more componets doesn't affect much, but small number of components seemed less stable. Number of components could be chosen based on elbow plot, like suggested here: https://cole-trapnell-lab.github.io/monocle-release/Paul_dataset_analysis_final.html

  * Cells are sampled to N=1,000 per cell group.

  * For this analysis we used for each sample tumor cells of that sample and cells from 3 normal duct cell populations that were pooled across samples.

2. ```Monocle_colorBySubtype.20230503.R``` -- script to color cells ordered along the built trajectory by cell type.

3. ```Monocle_corrTFs_withPT.20230503.R``` -- script to calculate correlation between pseudotime from the trajectories bult (using only the tumor cells and the respective cell-of-origin) and between the TF scores from chromVAR.

   * for this analysis we used sample-level objects only.



### Idea is to show, that tumor cells of different subtypes have different cell of origin:

    1. Lum <-- Lum mature.

    2. Basal <-- Lum progenitors.


### It can be concluded from the observation of what normal duct population is the closest to the tumor population of particular PAM50 subtype.