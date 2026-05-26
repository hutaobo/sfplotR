\name{plot_cophenetic_heatmap}
\alias{plot_cophenetic_heatmap}
\title{Plot a cophenetic distance heatmap}
\usage{
plot_cophenetic_heatmap(
  matrix,
  matrix_name = NULL,
  output_dir = NULL,
  output_filename = NULL,
  figsize = c(9, 9),
  cellwidth = 11,
  cmap = "RdBu",
  linewidths = 0.5,
  annot = TRUE,
  sample = "Sample",
  xlabel = NULL,
  ylabel = NULL,
  show_dendrogram = TRUE,
  return_xlabels = FALSE,
  parse_label = FALSE
)
}
\arguments{
\item{matrix}{Numeric matrix to plot.}
\item{matrix_name}{Optional matrix type label, such as \code{"row_coph"} or \code{"col_coph"}.}
\item{output_dir}{Directory where the PDF should be written.}
\item{output_filename}{Optional PDF filename.}
\item{figsize}{Numeric vector of PDF width and height in inches.}
\item{cellwidth}{Cell width and height passed to \code{pheatmap::pheatmap}.}
\item{cmap}{Color palette name. Currently \code{"RdBu"} uses RColorBrewer.}
\item{linewidths}{Retained for interface compatibility.}
\item{annot}{Retained for interface compatibility.}
\item{sample}{Sample label used in the plot title and default filename.}
\item{xlabel}{Optional x-axis label.}
\item{ylabel}{Optional y-axis label.}
\item{show_dendrogram}{Logical; whether to cluster rows and columns.}
\item{return_xlabels}{Logical; whether to return plotted x labels.}
\item{parse_label}{Logical; whether to parse labels through \pkg{latex2exp}.}
}
\value{
Writes a PDF heatmap. If \code{return_xlabels = TRUE}, returns plotted x labels;
otherwise returns \code{NULL} invisibly.
}
\description{
Creates a PDF heatmap for a cophenetic distance matrix using \pkg{pheatmap}.
}
