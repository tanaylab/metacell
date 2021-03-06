% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mc_apx_qc_cv.r
\name{mcell_wgtmc_pred_corr}
\alias{mcell_wgtmc_pred_corr}
\title{Compute predictive value of MC cover as correlation of weighted MC averages of sc umis}
\usage{
mcell_wgtmc_pred_corr(mc_id, mat_id, graph_id, do_log = F, focus_cells = NULL)
}
\arguments{
\item{mc_id}{id of metacells to analyze}

\item{mat_id}{umi matrix}

\item{do_log}{should log transofmration be performed}

\item{focus_cells}{should correltion be computed to a specific set of cells}
}
\description{
Compute predictive value of MC cover as correlation of weighted MC averages of sc umis
}
