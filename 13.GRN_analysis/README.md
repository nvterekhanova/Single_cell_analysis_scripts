# SCENIC analysis of Gene Regulatory Networks (GRNs).

## `Scripts` dir contains scripts for running SCENIC on cluster.

## `Post_scenic` dir contains scripts for running post-SCENIC analyses.


## ```Scripts``` contains scripts for running pySCENIC pipeline on the sampled object. We run pySCENIC 10 times (iterations), and then aggregate results across iterations, requiring that TF and its targets should be found in at least 8 iterations (80% of runs). See here for discussions from the authors and users related to itL

  * https://github.com/aertslab/pySCENIC/issues/169

  * https://pyscenic.readthedocs.io/en/latest/faq.html#how-can-i-prioritize-the-target-genes-within-a-particular-regulon (item 2).