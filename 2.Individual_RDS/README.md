# Changes to the atac_pipeline_v3.3.R and multiome_pipeline_v1.2.R:

## 2021-10-19:
small fixes because of issues with GenomeInfoDb (when adding annotation); not affecting results

## 2021-08-29, updates:
  
  * use signac v.1.3 for processing. 

  * use -neglog10(qval) as a score. 

  * remove peaks on chrY, and any peaks with "N" in their sequence. 

  * remove line "pairs=pairs[pairs$score_1>=pairs$score_2,]".

    + This will add small # of peaks from some peak clusters with peak number>=3.

    + This is a part of iterative removal procedure: previously we removed peaks that ov
erlapped any other peaks with greater signal (https://www.nature.com/articles/nmeth.4401
); without that line we will keep those if the overlapping peak with greater signal was overlapping another peak with even greater signal, and was removed first.

    + This is illustrated here: https://www.archrproject.com/bookdown/the-iterative-overlap-peak-merging-procedure.html.

    + Example: cluster of 3 overlapping peaks with scores: score_1>score_2>score_3. Previously we would keep peak_1 with score_1, and now we are keeping peak_1 and peak_3).
  

## 2021-07-04, update: modified to be compatible with multiome processing:

1. Also to be consistent with Cellranger-atac v.2/Cellranger-arc v.2. 

2. Signac v.1.2.1 and Seurat v.4.0.3.

3. Remove peaks on non-standart chromosomes, and in blacklist regions.


## Mofidified on 2020-08-31 using updates from Signac v.1.0.0, see here changes compared to Signac v.0.2.5: https://cran.r-project.org/web/packages/Signac/news/news.html.

   * References based on https://satijalab.org/signac/articles/pbmc_vignette.html.
