# Here are the scripts that need to perform mutation mapping to the cells. It uses mutation MAF obtained from bulk, and also 10Xmapping tool that identifies mutations present in single cell reads.


## To do mutation mapping, you need to follow these steps:

1. Select_muts_forMAF.20220714.R -- prepare MAF-file with mutations of interest. I just use all that were called in at least one sample of interest. Also, make sym-links to the .bam and .bam.bai.

2. getBarcodes.20220702.R -- extract barcodes-list for each sample, and save it in a folder. It will be used for creating mutation per cell matrix.

3. Map_muts.sh -- run 10Xmapping (Song's pipeline), step1/2 only.

4. create_matrix.20220702.pl -- create matrix using results from steps 2 & 3.

5. Extract_counts_forSeuratClusters.20220714.R -- aggregate results per cell group/cluster, and also intersect it with the WXS-calls, because we want to focus only on the mutations that were detected for this particular sample.