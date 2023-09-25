# CICERO analysis of chromatin co-accessibility:

  * Reference: https://www.sciencedirect.com/science/article/pii/S1097276518305471.


## Installation notes:

1. Need to install monocle3: https://cole-trapnell-lab.github.io/monocle3/docs/installation/.


2. To install it need to install spdep via conda: https://anaconda.org/conda-forge/r-spdep (only conda installtion worked for me)


3. Then you can go on and install CICERO:

```remotes::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")```
