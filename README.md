# nocurve

**Contextual empirical Bayes curve estimation**

The R package implements NOC, the proposed noise-orthogonal contextual estimator.
EB0 supplies the common empirical Bayes baseline. RAW uses unprojected
neighborhoods and provides the matched comparison for assessing projection.

This repository contains the R package and code for the manuscript's four
simulation scenarios and the single Raman training–test split.

## Installation and manuscript runs

Use R 4.5.2 for the simulation or R 4.6.0 for the Raman application. A C++17
toolchain and Python 3 are needed. Run the commands from the repository root.
The [dependency record](reproduce/dependency_pins.json) gives exact R package
versions. Install the pinned dependencies in an R library for the reproduction,
then install `nocurve` from this checkout:

```sh
R CMD INSTALL .
Rscript --vanilla reproduce/run.R simulation results-simulation
Rscript --vanilla reproduce/run.R raman results-raman
```

The simulation generates 100 independent noise draws per scenario using the
recorded seed rule and fits all eight methods to each input. It writes fitted
values, tuning parameters and original-grid loss summaries, including Monte
Carlo standard errors. The Raman script fits all eight methods to the training
spectrum and then reads the test reference to calculate RMSE, MAE and MAXAE.
Each command takes a new output directory.

The manuscript's simulation tables use the complete validated run on x86-64
Linux with R 4.5.2 and mgcv 1.9-4. This run covers 400 data sets and 3,200
method fits. For numerical comparisons with those tables, use the original
platform and the pinned dependency versions. Native REML fits can show
small numerical differences across platforms.

The [simulation design](reproduce/data/simulation.json) records the four mean
functions. The [Raman data note](reproduce/data/ATTRIBUTION.md) documents the
single split and attribution.
