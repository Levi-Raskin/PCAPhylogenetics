This is the GitHub repo for Raskin et al. (2025) Principal Components Analysis fails to recover phylogenetic structure in hominins.

This code works and should be sufficient to reproduce our analysis, but be forewarned that this code is not guaranteed in any way. Please reach out to [Levi Raskin](mailto:levi_raskin@berkeley.edu) if there are any issues and I would be very happy to help.

It has 5 subfolders.

data
  - contains the nexus file (from Mongle et al. 2023 10.1016/j.jhevol.2022.103311) used to infer the phylogenetic posterior distribution

scripts
  - DataAnalysis.R
    - Traditional morphological data simulations, PCA, and distance calculations
    - LDDMM simulated shapes read in, PCA, and distance calculations
    - Traditional morphological data analysis (mann whitney U tests, etc.)
    - LDDMM data analysis (mann whitney U tests, etc.)
    - Null distribution generation for the distance between two randomly generated trees and testing against the null
  - Figures.R
    - Code needed to reproduce the figures in the manuscript
  - ShapeSimulationCode
    - C++ files needed to perform the simulations of geometric morphometric data under LDDMM
  - mkUnordered.Rev
    - RevBayes file for posterior distribution generation

resultsGit
  - includes all the files needed to reproduce the results we got in our paper. Does not include raw simulated data or the raw tree file as that would be too large to store in github. Reach out to [Levi Raskin](mailto:levi_raskin@berkeley.edu) for access to the raw data.

manuscript
  - contains figures as svg, jpg, and adobe illustrator files
  - contains SOM tables s1 and s2

rSPR
  - contains 4 separate rspr executables (https://github.com/cwhidden/rspr) so that we can parallelize the SPR distance calculation. Probably there exist better ways to do than our implementation.
