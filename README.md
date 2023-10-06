# Treatment Effect Variable Importance Measures

This repo contains reproduction code for simulation and the applied example for the paper "Variable importance measures for heterogeneous causal effects"
by Oliver Hines, Karla Diaz-Ordaz, and Stijn Vansteelandt

## Simulation study

Simulation code is contained in `R/simulation`. 
To run simulations related to Data Generating Process (DGP) 1 and 2, run `experiment_1.R` from the repo root. This script takes around 20 hours on an M1 Macpro.
For DGP 3 run `experiment_2.R`. This takes around 6 hours on an M1 Macpro.
Simulation data, log files and plots will be saved to `Output/`.

## Illustrated example

The illustrated example uses AIDS Clinical Trials Group Protocol 175 (ACTG175) data which is obtained from the `speff2trial` package on CRAN. This can be run using `R/actg175/data_example.R`. The use of the super learner (and some imperfect parallelisation for Algorithm 2) means this example took around 30-60 mins to run on an M1 Macpro.

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
