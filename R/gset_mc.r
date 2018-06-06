#' adding genes to a gene set
#'
#' @param gset_id the basic gene set
#' @param genes names of genes to add
#' @param subset_id id of subset for adding the genes, 1 by defulat
#'
#' @export

mcell_gset_add_gene = function(gset_id, genes, subset_id = 1)
{
	gs = scdb_gset(gset_id)
	if(is.null(gs)) {
		stop("MC-ERR: trying to add gene to non existing gset id ", gset_id)
	}
	gs = gset_add_genes(scdb_gset(gset_id), genes, subset_id)
	scdb_add_gset(gset_id, gs)
}

#' Generate a gene set of markers from a metacell cover
#'
#' This should be mostly used for visualization and highlighting key genes
#'
#' @param gset_id id of gene set to generate
#' @param mc_id id of metacell object
#' @param filt_gset_id gene set to select from (null by default)
#' @param blacklist_gset_id gene set to exclude (null by default)
#'
#' @export
mcell_gset_from_mc_markers = function(gset_id, mc_id, 
					filt_gset_id=NULL, blacklist_gset_id = NULL)
{
	k_per_clust = get_param("scm_mc_mark_k_per_clust")
	min_gene_fold = get_param("scm_mc_mark_min_gene_fold")
	min_gene_cov = get_param("scm_mc_mark_min_gene_cov")
	
	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
		stop("MC-ERR: undefined metacell ", mc_id, " when selecting markers")
	}
	gene_folds = mc@mc_fp
	genes_pool = rownames(gene_folds)
	if(!is.null(filt_gset_id)) {
		filt_gset = scdb_gset(filt_gset_id)
		if(is.null(filt_gset)) {
			stop("MC-ERR: undefined filt gset id ",filt_gset_id)
		}
		genes_pool = intersect(genes_pool, names(filt_gset@gene_set))
		if(length(genes_pool) < 2) {
			stop("MC-ERR: intersect between metacell gene names and filt gene set is < 2")
		}
	}
	if(!is.null(blacklist_gset_id)) {
		bl_gset = scdb_gset(blacklist_gset_id)
		if(is.null(bl_gset)) {
			stop("MC-ERR: undefined blacklist gset id ",filt_gset_id)
		}
		genes_pool = setdiff(genes_pool, names(bl_gset@gene_set))
		if(length(genes_pool) < 2) {
			stop("MC-ERR: gene pool after blacklisting is < 2")
		}
	}

#	mask_folds = gene_folds[genes_pool,]*(mc@cov_gc[genes_pool,]>min_gene_cov)
	lfp = abs(log2(gene_folds[genes_pool,]))

	marks = unique(as.vector(unlist(
			apply(lfp,
				2,
				function(x)  {
				   names(head(sort(-x[x>min_gene_fold]),n=k_per_clust)) })
		     )))

	if(is.null(marks) | length(marks) < 2) {
		stop("MC-ERR no marks found, consider relaxing min_gene_cov or min_gene_fold")
	}

	smooth_n = max(3, round(ncol(lfp)/30))
	
	lfp_smoo = t(apply(lfp[marks,], 1, function(x) rollmean(x,smooth_n, fill=0)))
	lfp_smoo[,1:smooth_n] = rowMeans(lfp[marks,1:smooth_n])
	lfp_smoo[,ncol(lfp)-(1:smooth_n)] = 
				rowMeans(lfp[marks,ncol(lfp)-(1:smooth_n)])

	marks = marks[order(apply(lfp_smoo, 1, which.max))]

	marks_gs = rep(1, length(marks))
	names(marks_gs) = marks
	scdb_add_gset(gset_id, tgGeneSets(marks_gs, desc=sprintf("%s mc marks", mc_id)))
}
