\name{transcript_by_cell_analysis}
\alias{transcript_by_cell_analysis}
\title{Run transcript-by-cell spatial analysis}
\usage{
transcript_by_cell_analysis(
  cell_metadata,
  coords,
  sample_name = NULL,
  output_folder = NULL,
  coph_method = "average",
  n_jobs = 8
)
}
\arguments{
\item{cell_metadata}{Data frame containing at least \code{x}, \code{y}, and \code{celltype} columns.}
\item{coords}{Transcript coordinate data frame containing at least \code{x}, \code{y}, and \code{feature_name} columns.}
\item{sample_name}{Optional sample name used for output filenames.}
\item{output_folder}{Optional output directory.}
\item{coph_method}{Hierarchical clustering method passed to \code{stats::hclust}.}
\item{n_jobs}{Number of parallel workers.}
}
\value{
Returns \code{NULL} invisibly after writing Cell-GPS output files to disk.
}
\description{
Computes baseline cell-type topology and per-gene transcript-to-cell-type
cophenetic profiles, then writes heatmap and CSV outputs.
}
