#' Compute metacell using resampling iterations of graph cover k-means-like approach
#'
#' @param mc_id id of metacell object to be added
#' @param graph_id a knn graph object id into scdb
#' @param mat_id a matrix object id
#' @param min_mc_size target minimum metacell size. This is only an approximation and smaller MC may be returned by the algorithm
#' @param resamp_n number of resampling iterations
#'

mcell_coclust_from_graph_resamp = function(coc_id, 
		graph_id, 
		min_mc_size, p_resamp, n_resamp)
{
	tgs_clust_cool = get_param("scm_tgs_clust_cool")
	tgs_clust_burn = get_param("scm_tgs_clust_burn_in")

	graph = scdb_cgraph(graph_id)
	if(is.null(graph)) {
		stop("MC-ERR: cell graph id ", graph_id, " is missing when running add_mc_from_graph")
	}
	message("running bootstrap to generate cocluster")
	edges = graph@edges
	colnames(edges) = c("col1", "col2", "weight")

	K = round(nrow(edges)/length(graph@cell_names))

	resamp = tgs_graph_cover_resample(edges, knn = K, min_mc_size, cooling = tgs_clust_cool, burn_in = tgs_clust_burn, p_resamp = p_resamp, n_resamp = n_resamp)

	scdb_add_coclust(coc_id, 
			tgCoClust(graph_id=graph_id, coclust = resamp$co_cluster, n_samp=resamp$samples))
}