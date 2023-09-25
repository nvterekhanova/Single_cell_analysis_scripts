# Here are the scripts used for merging in pan-cancer snATAC project

* Select_peaks_for_pancan -- scripts for generating peaks used for merging. It identifies overlapping peaks, and prioritizes the ones of higher coverage. For this normalization of the peaks is also required, because coverage varies across samples.

* ../Merging_for_pancan_ATAC/Make_cohort_objects -- make cohort objects using pancancer set of peaks (so that they are the same across cancers).

* ../Merging_for_pancan_ATAC/Make_pancan_objects -- make pan-cancer merged objects using pan-cancer set of peaks.