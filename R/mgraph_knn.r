#' Compute metacell manifod graph using the confusion matrix of balanced K-nn between individual cells projected on metacells
#'
#' @param mgraph_id id of new object
#' @param mc_id meta cell id to work with
#' @param ignore_edges provide a data frame with mc1,mc2 pairs of edges to delete manually
#' @param graph_id graph_id of the similarity graph on cells from the metacell, in case confusion should be computed from similarities.
#' @param mctnetwork_id  id of a network object, in case confusion should be comuted from flows.
#' @param symetrize should the mc confusion matrix be symmetrized before computing layout?
#' @param ignore_mismatch try if the cgraph id can be partially overlapping with the metacell object - false by defualt and should be kept this way
#'
#' @export
mcell_mgraph_knn = function(mgraph_id, mc_id, 
			graph_id = NULL, mctnetwork_id = NULL,
			symmetrize=F, 
			ignore_mismatch = F,
			ignore_edges = NULL) 
{
	mc = scdb_mc(mc_id)
	if (is.null(mc)) {
		stop(sprintf("mc %s not found"), mc_id)
	}
	if(is.null(graph_id) & is.null(mctnetwork_id)) {
		stop("specfiy at least one of graph_id and mctnetwork_id when building an mgraph using Knn relations")
	}
	mgraph = mgraph_comp_knn(mc_id, graph_id, mctnetwork_id, 
					ignore_mismatch=ignore_mismatch, symmetrize=symmetrize)
	if(!is.null(ignore_edges)) {
		all_e = paste(mgraph$mc1, mgraph$mc2, sep="-")
		ig_e = paste(ignore_edges$mc1, ignore_edges$mc2, sep="-")
		ig_re = paste(ignore_edges$mc2, ignore_edges$mc1, sep="-")
		f = all_e %in% c(ig_e,ig_re)
		mgraph= mgraph[!f,]
		message("igoring ", sum(f), " edges")
	}
	scdb_add_mgraph(mgraph_id, tgMCManifGraph(mc_id, mgraph))
}

#' @export
mgraph_comp_knn = function(mc_id, graph_id, mctnetwork_id,
										ignore_mismatch=F, symmetrize=F)
{
	mgraph_K = get_param("mcell_mgraph_K")
	mgraph_T_edge = get_param("mcell_mgraph_T_edge")
	mgraph_max_confu_deg = get_param("mcell_mgraph_max_confu_deg")
	mgraph_edge_asym = get_param("mcell_mgraph_edge_asym")
	mgraph_k_expand_inout_factor = get_param("mcell_mgraph_expand_inout_factor")
	mgraph_max_fpcor_indeg = get_param("mcell_mgraph_max_fpcor_indeg")
	mgraph_max_fpcor_outdeg = get_param("mcell_mgraph_max_fpcor_outdeg")

	if(is.null(mgraph_max_confu_deg) & is.null(mgraph_max_fpcor_outdeg)) {
		stop("MC-ERR: Either max_confu_deg or max_fpcor_deg must be defined - currently both are null")
	}

	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
		stop("MC-ERR: mc id ", mc_id, " is missing when running add_mc_from_graph")
	}
	restrict_in_degree = T

	if(!is.null(graph_id)) {
		message("comp mc graph using the graph ", graph_id, " and K ", mgraph_K)
		confu = mcell_mc_confusion_mat(mc_id, graph_id, mgraph_K, 
													ignore_mismatch=ignore_mismatch)
	} else {
		message("comp mc graph using the flows ", mctnetwork_id, " and K ", mgraph_K)
		if(is.null(mctnetwork_id)) {
			stop("both graph id and network id are unspecified hen building knn mgraph")
		}
		mct = scdb_mctnetwork(mctnetwork_id)
		if(is.null(mct)) {
			stop("cannot get mctnetwork id ", mctnetwork_id, " in db")
		}
		confu = mctnetwork_get_flow_mat(mct, -1)
	}

	if(symmetrize) {
		confu =confu + t(confu)
	}
# k_expand_inout_factor=k_expand_inout_factor

	csize = as.matrix(table(mc@mc))
	csize = pmax(csize, 20)
	csize2 = csize %*% t(csize)
	csize2 = csize2 / median(csize)**2
	confu = confu / csize2

	confu_p_from = confu/rowSums(confu)
	confu_p_to = t(confu)/colSums(confu)

	if(!is.null(mgraph_max_confu_deg)) {
		rank_fr = t(apply(confu_p_from, 1, rank))
		rank_to = t(apply(confu_p_to, 1, rank))
		rank2 = rank_fr * rank_to
		diag(rank2) = 1e+6
		amgraph = apply(rank2, 1, function(x) {  rank(-x) <= (1+mgraph_max_confu_deg) })
		mgraph = amgraph * ((confu_p_from + confu_p_to) > mgraph_T_edge)
		if(restrict_in_degree) {
			amgraph2 = t(apply(rank2, 2, function(x) {  rank(-x) <= (1+mgraph_max_confu_deg) }))
			mgraph = mgraph * amgraph2
		}

		if(mgraph_edge_asym) {
			mgraph = amgraph * (confu_p_from>mgraph_T_edge)
			mgraph = amgraph * (t(confu_p_to)>mgraph_T_edge)
		}
		mgraph = mgraph>0 | t(mgraph>0)
	} else {
		mgraph = (confu_p_from + confu_p_to) > mgraph_T_edge
	}

	if(!is.null(mgraph_max_fpcor_outdeg)) {
		f = apply(abs(log2(mc@mc_fp)),1,max)>0.5
		fp_cor = tgs_cor(log2(mc@mc_fp[f,]))
		fp_rnk = t(apply(-fp_cor,1,function(x) rank(x)<=(mgraph_max_fpcor_outdeg+1)))
		fp_rnk_in = apply(-fp_cor,1,function(x) rank(x)<=(1+mgraph_max_fpcor_indeg))
		if(!is.null(mgraph_max_confu_deg)) {
				  mgraph = mgraph * fp_rnk * fp_rnk_in
		} else {
				  mgraph = fp_rnk * fp_rnk_in
		}
	}

	N = nrow(mgraph)
	e = which(mgraph>0)
	n1 = ceiling((e)/N)
	n2 = 1+((e-1) %% N)
	return(data.frame(mc1 = n1, mc2 = n2, dist=1))
}

