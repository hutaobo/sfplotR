# sfplotR

`sfplotR` is an R package for spatial structure analysis and visualization in spatial transcriptomics data. It provides R implementations of the SFplot / Cell-GPS workflow used to compute cluster-level cophenetic structure matrices, generate StructureMap-style heatmaps, and analyze transcript-to-cell spatial relationships.

This repository is maintained as a reviewer-friendly code companion for the associated manuscript.

## Main features

- Compute nearest-neighbor distance summaries between spatial clusters.
- Derive row-wise and column-wise cophenetic distance matrices from coordinate tables.
- Generate publication-style heatmaps for StructureMap and related outputs.
- Run transcript-by-cell analysis from cell metadata and transcript coordinates.
- Provide example scripts under `inst/example/` for manuscript-oriented analyses.
- Launch an interactive Shiny application for exploratory use.

## Repository layout

- `R/`: package source code.
- `man/`: package manual pages.
- `inst/example/`: example analysis scripts.
- `DESCRIPTION`: package metadata and dependencies.

## Installation

Install from GitHub with `remotes`:

```r
install.packages("remotes")
remotes::install_github("hutaobo/sfplotR")
```

Or install from a local clone:

```r
install.packages("devtools")
devtools::install("path/to/sfplotR")
```

The package targets R 4.3 or later.

## Minimal example

```r
library(sfplotR)

df <- data.frame(
  x = c(0, 1, 5, 6),
  y = c(0, 1, 5, 6),
  Cluster = c("A", "A", "B", "B")
)

result <- compute_cophenetic_distances_from_df(
  df = df,
  cluster_col = "Cluster",
  x_col = "x",
  y_col = "y"
)

plot_cophenetic_heatmap(
  matrix = result$row_cophenetic_df,
  matrix_name = "row_coph",
  sample = "Example",
  output_dir = "output"
)
```

## Core exported functions

- `compute_cluster_average_nn_distance_matrix`
- `compute_cluster_nn_distance_df`
- `compute_cophenetic_distances_from_df`
- `plot_cophenetic_heatmap`
- `transcript_by_cell_analysis`

## Notes for reviewers

- The main reviewer-facing code is in `R/`.
- Example scripts used during analysis development are available under `inst/example/`.
- Large local database artifacts are not part of the repository history and are intentionally ignored.
- A short repository walkthrough is available in [REVIEWER_GUIDE.md](REVIEWER_GUIDE.md).

## Citation

If you use this repository in connection with manuscript review, please cite the associated manuscript:

> Cophenetic Spatial Topology Embedding reveals multiscale tissue architecture in spatial omics

## License

This project is distributed under the MIT License.
