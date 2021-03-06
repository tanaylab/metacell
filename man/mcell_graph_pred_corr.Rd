% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mc_apx_qc.r
\name{mcell_graph_pred_corr}
\alias{mcell_graph_pred_corr}
\title{Compute predictive value of MC cover as correlation of MC averages of sc umis}
\usage{
mcell_graph_pred_corr(
  mat_id,
  graph_id,
  do_log = F,
  focus_cells = NULL,
  filt_genes = NULL
)
}
\arguments{
\item{mat_id}{id of matrix to analyze}

\item{graph_id}{id of cgraph to use for computing neighborhood per cell}

\item{do_log}{should log transofmration be performed}

\item{focus_cells}{should correltion be computed to a specific set of cells}
}
\description{
Compute predictive value of MC cover as correlation of MC averages of sc umis
}
