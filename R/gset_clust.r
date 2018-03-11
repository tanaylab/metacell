#' split gene set into several modules using clustering of genes by correlation over cells downsampled umi vector
#'
#' @param gset_id current gene set object. Only the list of genes will be used
#' @param mat_id scRNA matrix to use. We assume genes in the set are represented in the matrix, but if force is T we allow filtering of genes by represented matrix entry in addtion to clsutering
#' @param K how many gene modules to generate (kmeans parameter)
#' @param force if this is T, we allow clustering usign matrix in which not all genes are represnted. This will override the existing gset will a subset of the genes!
#' @param (ENV) scm_k_nonz_exp - to define transformation of UMIs
#'
#' @export

mcell_gset_split_by_dsmat = function(gset_id, mat_id, K, force=F)
{
	gset = scdb_gset(gset_id)
	if(is.null(gset)) {
		stop("missing gset id ", gset_id, " when trying to split gset into clusts")
	}	
	feat = gset_get_feat_mat(gset_id, mat_id, downsamp=T)

	gs  = names(gset@gene_set)	
	valid_gs = intersect(rownames(feat), gs)
	if(!force & length(valid_gs) != length(gs)) {
		stop("MC-ERR Missing genes in the matrix when trying to split gset by mat. Use force=T to filter these genes from the gene set before splitting genes into groups")
	}
	if(length(valid_gs) < 5) {
		stop("MC-ERR less than 5 genes when trying to split gset, this does not make sense, breaking")
	}
	k_nonz_exp = get_param("scm_k_nonz_exp")

	feat = log2(1+k_nonz_exp*as.matrix(feat))
	feat_cor = tgs_cor(t(feat))
	set.seed(19)
	hc = hclust(as.dist(1-feat_cor), "ward.D2")
	clst = cutree(hc, K)

	desc = sprintf("%s hc K=%d", gset@description,	K)
	scdb_add_gset(gset_id, gset_new_gset(clst, desc))
}

#' add new gene set based on an existing one with filtering specific clusters
#'
#' @param gset_id current gene set object. Only the list of genes will be used
#' @param filt_clusts list of clusters for filtering
#' @param new_id id of newly generated gene set
#'
#' @export

mcell_gset_remove_clusts = function(gset_id, filt_clusts, new_id, reverse=F)
{
	gset = scdb_gset(gset_id)
	if(is.null(gset)) {
		stop("trying to filter non existing gset, id ", gset_id)
	}
	sets = gset@gene_set
	if(!reverse) {
			  filt_sets = sets[!sets %in% filt_clusts]
	} else {
			  filt_sets = sets[sets %in% filt_clusts]
	}
	desc = sprintf("%s cl_filt", gset@description)
	scdb_add_gset(new_id, gset_new_gset(filt_sets, desc))
}
