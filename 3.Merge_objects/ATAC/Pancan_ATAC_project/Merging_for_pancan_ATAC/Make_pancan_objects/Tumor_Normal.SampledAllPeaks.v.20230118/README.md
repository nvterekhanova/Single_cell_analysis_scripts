# Here are the scripts for generating pan-cancer merged object that includes cancer cells and selected normal cell types (potential cells-of-origin).

1. ```Dim_reduction.CellsSampled.20221217.R``` -- script for merging, normalization and dimensional reduction.

2. ```RunChromVar.V2.R``` -- script to run chromVar to generate TF activity scores based on their motif accessibilites.


* ```archived/Filter_peaks.20230118.R``` -- script for peaks prioritization (was used in the previous mergings). Not used anymore, because this approach results in too small N of peaks.