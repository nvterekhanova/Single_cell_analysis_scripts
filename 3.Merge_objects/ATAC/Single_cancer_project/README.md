# Here are the scripts used for merging snATAC samples for BRCA project. Both snATAC and multiome samples were used.

## Steps:

1. ```Merge_ATAC_samples_prelim.v.1.0.R``` -- first, use this one to make initial merge using small # of peaks

2. ```Overlap_peaks_acrossSamples.20211016.R``` -- next, use this one to remove overlapping peaks, and requantify coverage in the new set of peaks.


  * Also, check scripts used for PanCan snATAC project: merging was done first on cohort level, and then on PanCancer level.

  * Peaks were also re-calculated for cohort objects first, and then for PanCan object.