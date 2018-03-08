#' metacell
#'
#' @import dplyr
#' @import ggplot2
#' @import data.table
#' @import tgconfig
#' @import cowplot
#' @import tgstat
#' @import tgutil
#' @import igraph
#' @importFrom pdist pdist
#' @importFrom cluster silhouette
#' @name metacell
#' @docType package

.onLoad <- function(libname, pkgname) {
	tgconfig::register_params(system.file('config/metacell_params.yaml', package='metacell'), package='metacell')
	ggplot2::theme_set(ggplot2::theme_classic())
}
