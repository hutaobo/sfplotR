\name{compute_cluster_average_nn_distance_matrix}
\alias{compute_cluster_average_nn_distance_matrix}
\title{Compute average nearest-neighbor distances between clusters}
\usage{
compute_cluster_average_nn_distance_matrix(
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
A data frame whose rows are source clusters and columns are target clusters,
with entries equal to mean nearest-neighbor distances.
}
\description{
Computes, for each source cluster, the average distance from cells in that
cluster to their nearest cell in each target cluster.
}
