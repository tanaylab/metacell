#' Build a cell graph using balanced knn after normalizing over a given metacell model
#'
#' @param mat_id matrix object id
#' @param gset_id gset object id defining the gene features used for computing distances
#' @param graph_id  new graph id to create
#' @param K the guideline Knn parameter. The balanced will be constructed aiming at K edges per cell
#'
#' @export

mcell_add_cgraph_bknn_norm_mc = function(graph_id, mat_id, gset_id, K, mc_to_norm)
{
	k_alpha = get_param("scm_balance_graph_k_alpha")
	k_beta = get_param("scm_balance_graph_k_beta")

	norm_mc = scdb_mc(mc_to_norm)
	if(is.null(norm_mc)) {
		stop("missing metacell for cgraph normalization, id = ", mc_to_norm)
	}
	feat = gset_get_feat_mat(gset_id, mat_id, downsamp=F)
	k_nonz_exp = get_param("scm_k_nonz_exp")
	feat = log2(1+k_nonz_exp*as.matrix(feat))

	raw_cor = tgs_cor(feat)

	raw_cor[is.na(raw_cor)] = 0

	N = ncol(feat)
	all_mc = rep(1,N)
	names(all_mc) = colnames(feat)
	all_mc[names(norm_mc@mc)] = norm_mc@mc+1

	n_out = N - length(norm_mc@mc)
	mc_size = tabulate(norm_mc@mc)
	w_cm = diag(c(1/n_out, 1/mc_size))[all_mc,]

	cor_mm = t(w_cm) %*% raw_cor %*% w_cm
	norm_cor = raw_cor - cor_mm[all_mc, all_mc]

	gr = cormat_to_balanced_knn(norm_cor, K, k_alpha, k_beta)

	cnames = colnames(feat)
	colnames(gr) = c("mc1","mc2","w")

	scdb_add_cgraph(graph_id, tgCellGraph(gr, cnames))
}

#' @export
cormat_to_balanced_knn = function(cormat, k_knn, k_alpha, k_beta)
{
	nodes = colnames(cormat)

	N = length(nodes)

	k_knn_potential = min(k_knn*k_alpha, N)

	m_knn = matrix(k_knn_potential, 
			nrow = N,ncol = N, 
			dimnames = list(nodes, nodes))
	message("fill m_knn with knn_ord")
	for(i in 1:N) {
		best_hits = head(order(cormat[i,]),k_knn_potential)
		m_knn[i,best_hits] = 1:k_knn_potential
	}

	message("start balancing ", N)
	M = k_knn*k_knn*10

	m_knn_io = pmax(-m_knn * t(m_knn) + M,0)

	A = nrow(m_knn_io)-k_knn*3	#no more than k_nn*3 contribs per row, or incoming neighbots
	m_knn_io_rows = t(apply(m_knn_io, 1, function(x) pmax(rank(x)-A,0)))
	message("done balance 1", N)
	A = nrow(m_knn_io)-k_knn	#no more than k_nn contribs per column, or outgoing neighbors
	m_knn_io = apply(m_knn_io_rows, 2, function(x) pmax(rank(x)-A,0))
	message("done balance 2", N)

	ij = as.data.frame(which(m_knn_io > 0, arr.ind=T))

	gr = data.frame(ij$row, ij$col, weight=m_knn_io[which(m_knn_io>0)])
	return(gr);
}
