% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/monocleHeatmap.R
\name{plot_pseudotime_heatmap2}
\alias{plot_pseudotime_heatmap2}
\title{Plots a pseudotime-ordered, row-centered heatmap which is slightly modified in monocle2}
\usage{
plot_pseudotime_heatmap2(
  cds_subset,
  cluster_rows = TRUE,
  hclust_method = "ward.D2",
  num_clusters = 6,
  hmcols = NULL,
  add_annotation_row = NULL,
  add_annotation_col = NULL,
  show_rownames = FALSE,
  use_gene_short_name = TRUE,
  norm_method = c("log", "vstExprs"),
  scale_max = 3,
  scale_min = -3,
  trend_formula = "~sm.ns(Pseudotime, df=3)",
  return_heatmap = FALSE,
  cores = 1
)
}
\arguments{
\item{cds_subset}{CellDataSet for the experiment (normally only the branching
genes detected with branchTest)}

\item{cluster_rows}{Whether to cluster the rows of the heatmap.}

\item{hclust_method}{The method used by pheatmap to perform hirearchical
clustering of the rows.}

\item{num_clusters}{Number of clusters for the heatmap of branch genes}

\item{hmcols}{The color scheme for drawing the heatmap.}

\item{add_annotation_row}{Additional annotations to show for each row in the
heatmap. Must be a dataframe with one row for each row in the fData table of
cds_subset, with matching IDs.}

\item{add_annotation_col}{Additional annotations to show for each column in
the heatmap. Must be a dataframe with one row for each cell in the pData table
of cds_subset, with matching IDs.}

\item{show_rownames}{Whether to show the names for each row in the table.}

\item{use_gene_short_name}{Whether to use the short names for each row. If
FALSE, uses row IDs from the fData table.}

\item{norm_method}{Determines how to transform expression values prior to
rendering}

\item{scale_max}{The maximum value (in standard deviations) to show in the
heatmap. Values larger than this are set to the max.}

\item{scale_min}{The minimum value (in standard deviations) to show in the
heatmap. Values smaller than this are set to the min.}

\item{trend_formula}{A formula string specifying the model used in fitting
the spline curve for each gene/feature.}

\item{return_heatmap}{Whether to return the pheatmap object to the user.}

\item{cores}{Number of cores to use when smoothing the expression curves shown
in the heatmap.}
}
\value{
A list of heatmap_matrix (expression matrix for the branch committment),
ph (pheatmap heatmap object),
annotation_row (annotation data.frame for the row), annotation_col (annotation
data.frame for the column).
}
\description{
The function plot_pseudotime_heatmap takes a CellDataSet object
(usually containing a only subset of significant genes) and generates smooth
expression curves much like plot_genes_in_pseudotime.
Then, it clusters these genes and plots them using the pheatmap package.
This allows you to visualize modules of genes that co-vary across pseudotime.
}
