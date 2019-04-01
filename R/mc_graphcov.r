#' Compute metacell using a native implementation of a graph cover k-means-like approach
#'
#' @param mc_id id of metacell object to be added
#' @param graph_id a knn graph object id into scdb
#' @param mat_id a matrix object id
#' @param min_mc_size target minimum metacell size. This is only an approximation and smaller MC may be returned by the algorithm
#'

mcell_add_mc_from_graph = function(mc_id, graph_id, mat_id, min_mc_size)
{
    old_seed = .set_seed(get_param("mc_rseed"))

	tgs_clust_cool = get_param("scm_tgs_clust_cool")
	tgs_clust_burn = get_param("scm_tgs_clust_burn_in")

	graph = scdb_cgraph(graph_id)
	if(is.null(graph)) {
		stop("MC-ERR: cell graph id ", graph_id, " is missing when running add_mc_from_graph")
	}
	mat = scdb_mat(mat_id)
	if(is.null(mat)) {
		stop("MC-ERR: mat id ", mat_id, " is missing when running add_mc_from_graph")
	}
	message("running graph clustering now - one iteration no bootstrap")
	edges = graph@edges
	colnames(edges) = c("col1", "col2", "weight")
	node_clust = tgs_graph_cover(edges, min_mc_size, cooling = tgs_clust_cool, burn_in = tgs_clust_burn)
	f_outlier = (node_clust$cluster == 0)

	outliers = colnames(mat@mat)[node_clust$node[f_outlier]]
	mc = as.integer(as.factor(node_clust$cluster[!f_outlier]))
	names(mc) = colnames(mat@mat)[!f_outlier]
	message("building metacell object, #mc ", max(mc))
	cell_names = colnames(mat@mat)
	scdb_add_mc(mc_id, tgMCCov(mc, outliers, mat))
	message("reordering metacells by hclust and most variable two markers")

	.restore_seed(old_seed)

	mcell_mc_reorder_hc(mc_id)
}
