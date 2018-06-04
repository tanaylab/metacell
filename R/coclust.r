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

#' Return a filter (boolean vector) selecting only coclust edges that are nearly as frequent as a user defined K-nn parameter
#' 
#' @param coc_id coclust object id
#' @param K the number of top-coclustering neighbors to consider per node
#' @param alpha the relexation paramter to apply for filtering coclustering neighbors below the ttop K ones. A pair (n1,n2) with weight w will be filtered if knn(n1,K) > w*alpha or knn(n2,K) > w*alpha
#'
#' @export

mcell_coclust_filt_by_k_deg = function(coc_id, K, alpha)
{
	coc = scdb_coclust(coc_id)
	if(is.null(coc)) {
		stop("MC-ERR: coclust ", coc_id , " is missing when trying to derive K-deg filter")
	}

	edges = coc@coclust
	
	deg_wgt = as.matrix(table(c(edges$node1, edges$node2), c(edges$cnt,edges$cnt)))
	deg_cum = t(apply(deg_wgt, 1, function(x) cumsum(rev(x))))
	thresh_Kr = rowSums(deg_cum > K)
	thresh_K = rep(NA, length(levels(edges$node1)))
	names(thresh_K) = levels(edges$node1)
	thresh_K[as.numeric(names(thresh_Kr))] = thresh_Kr

	filt_edges = thresh_K[edges$node1] < edges$cnt * alpha | 
							thresh_K[edges$node2] < edges$cnt * alpha

	return(filt_edges)
}

