#' manifold graph structure over a metacell object
#'
#' Splitting metacells over a discrete time axis, defining manifold connections and estimated flows over them
#'
#' @slot mc_id id of the metacell object we represent as a network
#' @slot times_nms names of the time points (Default 1:T)
#' @slot mc_t distribution of metacells (rows) over time points (cols)
#' @slot mc_manifold a data frame defining triplets mc1, mc2, distance.
#'
#' @export tgMCManifGraph
#' @exportClass tgMCManifGraph
tgMCManifGraph <- setClass(
   "tgMCManifGraph",
	slots = c(
	  mc_id = "character",
	  mgraph = "data.frame"
	)
)

#' Construct a meta cell manifold graph
#'
#'
#' @param mc_id metacell object id
#' @param mgraph data fra,e defining mc1, mc2, distance
#' @export

setMethod(
  "initialize",
  signature = "tgMCManifGraph",
  definition =
    function(.Object, mc_id, mgraph) {
		.Object@mc_id = mc_id
		.Object@mgraph = mgraph 
		mc = scdb_mc(mc_id)
		if(is.null(mc)) {
			stop("MC-ERR unkown mc_id ", mc_id, " when building mc mgraph")
		}
      return(.Object)
    }
)

#' Generate a new metacell manifold graph object
#'
#' This constructs a meta cell manifold graph object - only encapsulating an edge list data frame
#'
#' @param mc_id id of scdb meta cell object ot be added
#' @param mgraph the mgraph data frame containing fields mc1, mc2, distance
#' @export
mcell_new_mc_mgraph = function(mgraph_id, mc_id, mgraph)
{
	scdb_add_mc(mgraph_id, tgMCManifGraph(mc_id, mgraph))
}
