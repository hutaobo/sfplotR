\name{compute_cluster_nn_distance_df}
\alias{compute_cluster_nn_distance_df}
\title{Compute cell-to-cluster nearest-neighbor distances}
\usage{
compute_cluster_nn_distance_df(
  df,
  cluster_col = "Cluster",
  x_col = "x",
  y_col = "y"
)
}
\arguments{
\item{df}{A data frame containing cluster labels and spatial coordinates.}
\item{cluster_col}{Name of the column containing cluster labels.}
\item{x_col}{Name of the x-coordinate column.}
\item{y_col}{Name of the y-coordinate column.}
}
\value{
A data frame with one row per cell, one column per target cluster, and a
\code{cell_cluster} column containing each cell's source cluster.
}
\description{
Computes each cell's nearest-neighbor distance to every cluster in the input
data frame.
}
