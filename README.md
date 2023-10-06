# Treatment Effect Variable Importance Measures

This repo contains reproduction code for simulation and the applied example for the paper "Variable importance measures for heterogeneous causal effects"
by Oliver Hines, Karla Diaz-Ordaz, and Stijn Vansteelandt

## Simulation study

Simulation code is contained in `R/simulation`. To run simulations, run `experiment_0.R` from the repo root. Simulation data, log files and plots will be saved to `Output/`.

## Illustrated example

Todo: write readme.

## Dependencies

This code has been tested using `R v4.2.1`.
Please make sure the following CRAN dependencies are installed:

- `tidyverse`
- `arrow`
- `ranger`
- `mgcv`
- `glmnet`
- `gam`
- `xgboost`
- `lightgbm`
- `SuperLeaner`
- `cowplot`
- `latex2exp`
- `speff2trial`
