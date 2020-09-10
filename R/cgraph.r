#' Cell graph
#'
#' Representing a cell-to-cell similarity graph.
#'
#' @slot cell_names names of cells (Ids in edges relates to these)
#' @slot edges data frame specifiying the graph weighted edges as returned from tgs_cor_graph
#' @slot nodes vector of node names
#'
#' @export tgCellGraph
#' @exportClass tgCellGraph
tgCellGraph <- setClass(
   "tgCellGraph",
	slots = c(
	  cell_names = "vector",
	  edges	= "data.frame",
	  nodes = "vector")
)

#'
#'
#' This constructs graph object from a data frame on columns col1,col2, weight, (the format returened by tgs_cor_graph). Note this is not intended to serve natively complex graph algs etc, but is just a container for use by scdb and operaitons on scmats (creating a graph) and mcell (using the graph to compute graph covers)
#'
#' @param gr data frame with columns col1 (first node), col2 (second node), weight - float between 0 and 1.
#' @export

setMethod(
  "initialize",
  signature = "tgCellGraph",
  definition =
    function(.Object, gr, cnames) {
		.Object@cell_names = cnames
		.Object@edges = gr
		n1 = unique(as.character(gr$mc1))
		n2 = unique(as.character(gr$mc2))
		.Object@nodes = unique(c(n1,n2))
		if(length(cnames) != length(.Object@nodes)) {
			miss = setdiff(cnames, .Object@nodes)
			message("sim graph is missing ", length(miss), " nodes, out of ", length(cnames))
		}
      return(.Object)
    }
)

