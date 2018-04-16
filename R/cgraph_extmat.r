#' Build a cell graph using blanacing of an extrenal distance matrix
#'
#' @param mat_id matrix object id
#' @param d_mat distance matrix, col and row names should be a super set of the non ignored cells for the matrix in mat_id
#' @param graph_id  new graph id to create
#' @param K the guideline Knn parameter. The balanced will be constructed aiming at K edges per cell
#' @param k_expand determine how much shoudl the K be expanded intially in order to find enough balanced neighbors.
#'
#' @export

mcell_add_cgraph_from_distmat= function(mat_id, d_mat, graph_id, K, k_expand=10)
{
	mat = scdb_mat(mat_id)
	if(is.null(mat)) {
		stop("cannot get matrix ", mat_id, " whne building a cgraph from d_mat")
	}
	cnames = intersect(colnames(mat@mat), colnames(d_mat))
	if(length(cnames) != ncol(mat@mat)) {
		stop("incomaptible cell (col names) in matrix ", mat_id, " and supplied distance matrix, missing ", length(setdiff(colnames(mat@mat), colnames(d_mat))), " cells in d_mat")
	}

	message("willl build balanced knn graph on ", ncol(d_mat), " cells, this can be a bit heavy for >20,000 cells")
	k_beta = get_param("scm_balance_graph_k_beta")

	d_mat = as.matrix(d_mat)
	x_knn = tgs_knn(max(d_mat)-d_mat, K*k_expand)
	gr = tgs_graph(x_knn, K, k_expand, k_beta)

	colnames(gr) = c("mc1","mc2","w")

	scdb_add_cgraph(graph_id, tgCellGraph(gr, cnames))
}
