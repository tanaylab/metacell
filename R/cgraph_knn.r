#' Build a cell graph using balanced knn graph on given gene features
#'
#' @param mat_id matrix object id
#' @param gset_id gset object id defining the features used for distances
#' @param graph_id  new graph id to create
#' @param K the guideline Knn parameter. The balanced will be constructed aiming at K edges per cell
#'
#' @export

mcell_add_cgraph_from_mat_bknn = function(mat_id, gset_id, graph_id, K, dsamp=F)
{
	feat = gset_get_feat_mat(gset_id, mat_id, downsamp=dsamp, add_non_dsamp=dsamp)
	k_nonz_exp = get_param("scm_k_nonz_exp")
	feat = log2(1+k_nonz_exp*as.matrix(feat))

	message("willl build balanced knn graph on ", ncol(feat), " cells and ", nrow(feat), " genes, this can be a bit heavy for >20,000 cells")
	k_alpha = get_param("scm_balance_graph_k_alpha")
	k_beta = get_param("scm_balance_graph_k_beta")
	gr = tgs_cor_graph(x=feat, knn=K, k_expand=10, k_alpha=k_alpha, k_beta=k_beta)
	cnames = colnames(feat)
	colnames(gr) = c("mc1","mc2","w")

	scdb_add_cgraph(graph_id, tgCellGraph(gr, cnames))
}
#' Compute a cell cell correlation matrix using features defined by a gene set
#'
#' @param mat_id matrix object id
#' @param gset_id gset object id defining the features used for distances
#'
#' @export

mcell_gen_cell_cor_gset = function(mat_id, gset_id, dsamp=F)
{
	feat = gset_get_feat_mat(gset_id, mat_id, downsamp=dsamp, add_non_dsamp=dsamp)
	k_nonz_exp = get_param("scm_k_nonz_exp")
	feat = log2(1+k_nonz_exp*as.matrix(feat))

	return(tgs_cor(feat))
}
