# Scripts for making cohort level objects:

1. ```snRNA```.

2. ```snATAC``` -- Requires making set of non-overlapping peaks across samples before merging (peaks harmonization).

  * Peak harmonization is required because we need to merge assays from different objects, and unlike snRNA objects that all have same rows corresponding to genes, snATAC objects will have different rows.

  * This is because peaks were called independantly for each object, and so they have different coordinates.