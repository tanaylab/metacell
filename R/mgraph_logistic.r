#' Compute metacell manifold graph using logistic distances and balanced K-nn
#'
#' @param mgraph_id id of new object
#' @param mc_id meta cell id to work with
#' @param ignore_edges provide a data frame with mc1,mc2 pairs of edges to delete manually
#' @param feats_gset gene set name for use for computing distances
#' @param feats_exclude list of genes to exclude from the features gene set
#' @param logist_loc the "location" parametr of the logistic function used to determine parametric distances between metacelles
#' @param logist_scale the "lscale" parametr of the logistic function used to determine parametric distances between metacelles
#'
#' @export
mcell_mgraph_logistic = function(mgraph_id, mc_id, 
			feats_gset,
			ignore_edges = NULL, 
			feats_exclude =NULL,
			logist_loc = 1, logist_scale = 0.2, 
			logist_eps = 4e-5, max_d_fold = 3)
{
	mc = scdb_mc(mc_id)
	if (is.null(mc)) {
		stop(sprintf("mc %s not found"), mc_id)
	}
	feat_genes = NULL
	gset = scdb_gset(feats_gset)
	if(is.null(gset)) {
		stop("Unkown gset ", feats_gset, " defined in mgraph logist construction")
	}
	feat_genes = names(gset@gene_set)
	if(!is.null(feat_genes)) {
		feat_genes = setdiff(feat_genes, feats_exclude)
	}
	message("got ", length(feat_genes), " feat genes for mc graph construction")

	mgraph = mgraph_comp_logist(mc, feat_genes, 
					logist_loc, logist_scale, logist_eps, max_d_fold)

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
mgraph_comp_logist= function(mc, genes, loc, scale, eps, max_d_fold)
{
	max_deg = get_param("mcell_mgraph_max_confu_deg")
	legc = log2(mc@e_gc[genes,] + eps)

	logist_d = function(x) {
		d = abs(legc - x)
		d = plogis(d, loc, scale)
		return(colSums(d))
	}
	a = apply(legc, 2, logist_d) 

#connect - d-best outgoing. filter by d_best_ratio < 2
	diag(a) = 1000;
	d_T = apply(a, 1, function(x) sort(x,partial=2)[2])
	a_n = a/d_T
	diag(a) = 0;
	diag(a_n) = 0;

   rank_fr = t(apply(a, 1, rank))
   rank_fr_m = rank_fr
   rank_fr_m[a_n > max_d_fold] = 1000
   rank_fr_m2 = rank_fr_m
   rank_fr_m2[t(a_n) > max_d_fold] = 1000

   edges = apply(rank_fr_m2, 1, function(x) {
                        mc2 = which(x > 0 & x <= max_deg+1);
                        mc1 = rep(which.min(x), length(mc2));
                        return(data.frame(mc1 = mc1, mc2=mc2)) })
   ed_df = as.data.frame(do.call(rbind, edges))
	ed_df$dist = apply(ed_df, 1, function(x) 1+a[x[1],x[2]])
	return(ed_df)
}
