# Reviewer Guide

This repository contains the R package companion to the `sfplot` / Cell-GPS workflow described in the manuscript.

## Suggested reading order

1. `README.md`
2. `R/compute_cophenetic_distances_from_df.R`
3. `R/plot_cophenetic_heatmap.R`
4. `R/tbc_analysis.R`
5. `inst/example/`

## What each area contains

- `R/`
  Core package implementation.
- `inst/example/`
  Example analysis scripts used during workflow development.
- `DESCRIPTION`
  Package metadata and declared dependencies.

## Minimal install

```r
install.packages("remotes")
remotes::install_github("hutaobo/sfplotR")
```

## Minimal input contract

For the basic cophenetic-distance workflow, the minimal required columns are:

- `x`
- `y`
- a cluster label column such as `Cluster`

## Reproducibility note

Large local GEO database files and external datasets are not committed to this repository. The repository is intended to expose the analysis code and package structure used in the study without bundling large local artifacts.
