#' Linear simple pipeline for turning a matrix to a MC and plotting std figs
#'
#' The simplification here is that all objects (mat, gstat, marker gset, cgraph, coclust, mc, mc2d) will have the same ID and will be produced linerly without much questions asked. Also the parameters are simplified for typical use in small datasets.
#'
#'	@param mat_id the id to start from
#' @param T_vm threshold on normalized varmin, for selecting markers. Lower this is you don't get enough info from the selected markers. (0.2)
#' @param T_tot threshold on minimal total expression for markers (200)
#' @param T_top3 threshold on minimal number of genes in the top3 expressing cells for markers (3)
#' @param Knn the raw Knn parameter used for constructing a balanced Knn graph. This should be set to roughly the size of desired metacells (120)
#' @param n_resamp number of bootstraps (500). Decrease to save time.
#' @param T_weight number of co-clustering occurences in bootstrap to add an edge in the phase of generating metacells from the bootstrap co-clustering matrix.
#'
#' @export

mcell_pipe_mat2mc2d = function(mat_id,
	T_vm=0.2,
	T_tot=200,
	T_top3=3,
	Knn=120,
	n_resamp=500,
	T_weight=round(n_resamp*0.75/8))
{
	id = mat_id
	id_plot_marks = sprintf("%s.plot_marks", id)

	mcell_add_gene_stat(id, id)

	mcell_gset_filter_varmean(id, id, T_vm=T_vm, force_new=T)
	mcell_gset_filter_cov(id, id, T_tot=T_tot, T_top3=T_top3)

	mcell_add_cgraph_from_mat_bknn(id, id, id, K=Knn)

	mcell_coclust_from_graph_resamp(id, id,
	                                min_mc_size=round(Knn/6),
	                                p_resamp=0.75,
	                                n_resamp=n_resamp)

	mcell_mc_from_coclust_louv_sub(id, id, id,
			 max_clust_size=4000,
			 max_mc_size=100,
			 min_mc_size=round(Knn/6),
			 T_weight = T_weight)

	mcell_gset_from_mc_markers(gset_id=id_plot_marks, mc_id=id)

	if(!is.null(mark_fn)) {
	  marks_colors = read.table(mark_fn, sep="\t", h=T, stringsAsFactors=F)
	  mc_colorize(mc_id=id, marker_colors=marks_colors)
	}

	mcell_mc_plot_marks(mc_id=id,
						gset_id=id,
						mat_id=id)

	mcell_mc2d_force_knn(id, id, id)
	mcell_mc2d_plot(mc2d_id=id)
}
