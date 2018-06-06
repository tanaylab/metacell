#' Metacell colustering object
#'
#' Represent coclustering data derived by resmapling iteration and graph cover (or any other coverage/clustering heuristic)
#'
#' @slot graph_id id of the scdb similairty graph object
#' @slot coclust tidy format matrix with co-occurence stats
#' @slot n_samp number of times each node was sampled
#'
#'
#' @export tgCoClust
#' @exportClass tgCoClust
tgCoClust <- setClass(
   "tgCoClust",
	slots = c(
	  graph_id = "character",
	  coclust = "data.frame",
	  n_samp = "vector")
)

#' Construct a coclust object
#' 
#'
#' @param graph_id just cahce the grpah id
#' @param coclust data frame with fields node1, node2, count
#' @param n_samp vector with statistics on number of time each node was sampled
#' @export

setMethod(
  "initialize",
  signature = "tgCoClust",
  definition =
    function(.Object, graph_id, coclust, n_samp) {
		if(is.null(scdb_cgraph(graph_id))) {
			stop("initlaizing coclust using graph id ", graph_id, " that is missing from scdb")
		}
		if(class(coclust)[1] != "data.frame") {
			stop("initializing coclust with non data frame coclust parameter")
		}
		cc_nms = colnames(coclust)
		if(!"node1" %in% cc_nms | !"node2" %in% cc_nms | !"cnt" %in% cc_nms) {
			stop("bad fields in coclust parameter")
		}
		.Object@graph_id= graph_id
		.Object@coclust = coclust
		.Object@n_samp = n_samp
      return(.Object)
    }
)
