#' gnereate/filter gene features from gstat normalized var/mean 
#'
#' @param gstat_id the ID of the gstat object to use
#' @param gset_id if this exists, filter_varmean will restrict hte current genes in the set with genes matching the var mean threshold, if not, it will generate a new gene sets object with one set including all high variance genes
#' @param T_vm the threshold on normalized var/mean, recommended values are usually around 0.2, but this may vary with the data
#'
#' @export

mcell_gset_filter_varmean = function(gstat_id, gset_id, T_vm, force_new=F)
{
	gstat = scdb_gstat(gstat_id)
	if(is.null(gstat)) {
		stop("missing gstat with id ", gstat_id, " when trying to generate a varmean gene set")
	}
	vm_genes = rownames(gstat)[gstat$ds_vm_norm > T_vm]
	gset = scdb_gset(gset_id)
	if(!is.null(gset) & force_new) {
		scdb_del_gset(gset_id)
		gset = NULL
	}
	if(is.null(gset)) {
		vm_set = rep(1, length(vm_genes))
		names(vm_set) = vm_genes
		gset = gset_new_gset(vm_set, sprintf("VM %f",T_vm))
		scdb_add_gset(gset_id, gset)
	} else {
		gset_vm = gset_new_restrict_gset(gset, vm_genes, desc=sprintf("%s VM %f", gset@description, T_vm))
		scdb_add_gset(gset_id, gset_vm)
	}
}

#' gnereate/filter gene features from coverage threshold in gstat table
#'
#' @param gstat_id the ID of the gstat object to use
#' @param gset_id if this exists, filter_varmean will restrict hte current genes in the set with genes matching the var mean threshold, if not, it will generate a new gene sets object with one set including all high variance genes
#' @param T_tot total down sampled coverage threhsold
#' @param T_top3 threshold value for the third highest umi count for the gene (only genes with top3>T_top3 are used)
#'
#' @export

mcell_gset_filter_cov = function(gstat_id, gset_id, T_tot, T_top3, force_new=F)
{
	gstat = scdb_gstat(gstat_id)
	if(is.null(gstat)) {
		stop("missing gstat with id ", gstat_id, " when trying to generate a varmean gene set")
	}
	cov_genes = rownames(gstat)[gstat$tot > T_tot & gstat$ds_top3 > T_top3]
	gset = scdb_gset(gset_id)
	if(!is.null(gset) & force_new) {
		scdb_del_gset(gset_id)
		gset = NULL
	}
	if(is.null(gset)) {
		cov_set = rep(1, length(cov_genes))
		names(cov_set) = cov_genes
		gset = gset_new_gset(cov_set, sprintf("Tot %d top3 %d",T_tot, T_top3))
		scdb_add_gset(gset_id, gset)
	} else {
		gset_cov = gset_new_restrict_nms(gset, cov_genes, desc=sprintf("%s tot %d top3 %d", gset@description, T_tot, T_top3))
		scdb_add_gset(gset_id, gset_cov)
	}
}

#' gnereate/filter gene features from statistics on correlation with umi count
#'
#' @param gstat_id the ID of the gstat object to use
#' @param gset_id if this exists, filter_varmean will restrict hte current genes in the set with genes matching the var mean threshold, if not, it will generate a new gene sets object with one set including all high variance genes
#' @param T_szcor upper limit on normalized sz_cor (low values mark interesting gene features). If you use this, consider values around -0.1 - but evaluate carefully your decision using the gstat empirical data
#'
#' @export

mcell_gset_filter_szcor = function(gstat_id, gset_id, T_szcor, force_new=F)
{
	gstat = scdb_gstat(gstat_id)
	if(is.null(gstat)) {
		stop("missing gstat with id ", gstat_id, " when trying to generate a varmean gene set")
	}
	szcor_genes = rownames(gstat)[gstat$sz_cor_norm < T_szcor]
	gset = scdb_gset(gset_id)
	if(!is.null(gset) & force_new) {
		scdb_del_gset(gset_id)
		gset = NULL
	}
	if(is.null(gset)) {
		szcor_set = rep(1, length(szcor_genes))
		names(szcor_set) = szcor_genes
		gset = gset_new_gset(szcor_set, sprintf("szcor %f",T_szcor))
		scdb_add_gset(gset_id, gset)
	} else {
		gset_szcor = gset_new_restrict_nms(gset, szcor_genes, desc=sprintf("%s szcor %f", gset@description, T_szcor))
		scdb_add_gset(gset_id, gset_cov)
	}
}

#' Select/filter gene features from using multiple statistics from the gstat table. All genes passing the selected thresholds are included
#'
#' @param gstat_id the ID of the gstat object to use
#' @param gset_id if this exists, the function will restrict the current genes in the set with genes matching the selected thresholds, if not, it will generate a new gene sets object with one set including all selected genes
#' @param T_tot total down sampled coverage thresholds (genes with tot UMIs < T_tot are filtered out)
#' @param T_top3 threshold value for the third highest umi count for the gene (genes with top3<T_top3 are filtered out)
#' @param T_szcor threshold value for the normalized size correlation statistic (only genes with sz_cor < T_szcor are selected). If you use this, consider values around -0.1 - but evaluate carefully your decision using the gstat empirical data
#' @param T_vm the threshold value for the normalized var/mean (only genes with varmean > T_vm are selected) Recommended values are usually around 0.2, but this may vary with the data. Not recommended for datasets with hihgly heterogeneous cell sizes (e.g. in whole-organisms datasets)
#' @param T_niche threshold value for the normalized niche score statistic (only genes with niche_norm > T_niche are selected). Recommended to use in combination with szcor to add genes with strongly restricted expression patterns. Consider using values around 0.05
#' @param blacklist option list of gene IDs to be excluded
#' @param force_new will overwrite existing gene set object (gset_id) in the database if it exists
#' @export

mcell_gset_filter_multi = function(gstat_id, gset_id, T_tot, T_top3,T_szcor=NULL, T_vm=NULL,T_niche=NULL, force_new=F,blacklist=c())
{
	gstat = scdb_gstat(gstat_id)
	if(is.null(gstat)) {
		stop("missing gstat with id ", gstat_id, " when trying to generate a varmean gene set")
	}
	
	if(!is.null(T_szcor)){szcor_genes = rownames(gstat)[gstat$sz_cor_norm < T_szcor]}else{szcor_genes=c()}
	if(!is.null(T_vm)){vm_genes = rownames(gstat)[gstat$ds_vm_norm > T_vm]}else{vm_genes=c()}
	if(!is.null(T_niche)){niche_genes = rownames(gstat)[gstat$niche_norm > T_niche]}else{niche_genes=c()}
	
	union_genes=union(szcor_genes,vm_genes)
	union_genes=union(union_genes,niche_genes)
	
	cov_genes = rownames(gstat)[gstat$tot > T_tot & gstat$ds_top3 > T_top3]
	filtered_genes=intersect(cov_genes,union_genes)
	
	filtered_genes=setdiff(filtered_genes,blacklist)
	message("Selected ",length(filtered_genes)," markers")
	
	gset = scdb_gset(gset_id)
	if(!is.null(gset) & force_new) {
		scdb_del_gset(gset_id)
		gset = NULL
	}

	gene_set = rep(1, length(filtered_genes))
	names(gene_set) = filtered_genes
	gset = gset_new_gset(gene_set, sprintf("szcor %f, vm %f, niche %f, tot %f, top3 %f",T_szcor,T_vm,T_niche, T_tot, T_top3))
	scdb_add_gset(gset_id, gset)
}