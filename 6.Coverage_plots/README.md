# Here are scripts for making coverage plots using snATAC/multiome

  * and also its possible to add expression track to the coverage plots using separate merged RNA object.

```Coverage_plots.SelectedEpigenetcallyDriven.20230329.R``` -- to generate coverage plots for pan-cancer snATAC paper (Extended Data Fig. 3c).

```Coverage_plots_PAM50_markers_Fig2E.R``` -- to generate coverage plots for HTAN BRCA paper (current Fig. 2E).



## Notes:

  * Use CoveragePlot and CombineTracks.

1. https://satijalab.org/signac/reference/coverageplot -- need to have coverage, expression and other tracks separately. They will be combined using ```CombineTracks```.

2. https://satijalab.org/signac/reference/combinetracks -- provide expression data here to the ```expression.plot``` arg.
