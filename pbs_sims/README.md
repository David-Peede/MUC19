# `pbs_sims`

This directory contains the code to replicate the PBS simulations.

## Packages

All packages are publicly available and their documentation can be viewed at the following places:

- [`SLiM v3.7.1`](https://messerlab.org/slim/)
- [`NumPy`](https://numpy.org/doc/stable/reference/index.html)
- [`pandas`](https://pandas.pydata.org/docs/)
- [`ggplot2`](https://ggplot2.tidyverse.org/)
- [`ggExtra`](https://cran.r-project.org/web/packages/ggExtra/vignettes/ggExtra.html)
- [`ggpubr`](https://rpkgs.datanovia.com/ggpubr/)

## Overview

- `data`
  - _MUC19_ annotation inputs for `SLiM`.
- `output`
  - Output files from `SLiM` where `snegative` represents simulations using the demographic model described in [Medina-Muñoz et. al., 2023](https://doi.org/10.1101/2023.03.06.531060) and `vnegative` represents simulations using the demographic model described in Figure S21.
- `scripts`
  - The code used to run the simulations.
  - A walkthrough of the pipeline can be found at `./scripts/sim_pipeline.ipynb`
