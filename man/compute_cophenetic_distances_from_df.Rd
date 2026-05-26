\name{compute_cophenetic_distances_from_df}
\alias{compute_cophenetic_distances_from_df}
\title{Compute cophenetic distance matrices from spatial metadata}
\usage{
compute_cophenetic_distances_from_df(
  df,
  cluster_col = "Cluster",
  x_col = "x",
  y_col = "y",
  method = "average",
  print = FALSE
)
}
\arguments{
\item{df}{A data frame containing cluster labels and spatial coordinates.}
\item{cluster_col}{Name of the column containing cluster labels.}
\item{x_col}{Name of the x-coordinate column.}
\item{y_col}{Name of the y-coordinate column.}
\item{method}{Hierarchical clustering method passed to \code{stats::hclust}.}
\item{print}{Logical; whether to print distance range diagnostics.}
}
\value{
A list with normalized row and column cophenetic distance matrices:
\code{row_cophenetic_df} and \code{col_cophenetic_df}.
}
\description{
Builds a cluster-by-cluster nearest-neighbor distance representation and
returns normalized cophenetic distance matrices from row and column
hierarchical clustering.
}
