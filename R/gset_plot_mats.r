#' Plot gene set correlation matrices given a an scamt. See version using metacells for potentially more robust behavior. This is used to detemrine initial feature selectio (e.g. filtering biologically irrelevant gene modules)
#'
#'
#'
#' @export

mcell_plot_gset_cor_mats = function(gset_id, scmat_id, downsamp=T)
{
	feat = gset_get_feat_mat(gset_id, scmat_id, downsamp=downsamp)
	gset = scdb_gset(gset_id)
	k_nonz_exp = get_param("scm_k_nonz_exp")
	feat = as.matrix(log2(1+k_nonz_exp*feat))

	mat = scdb_mat(scmat_id)
	cnms = colnames(feat)
	csize = colSums(mat@mat[,cnms])

	feat = rbind(feat, csize)
	rownames(feat)[nrow(feat)] = "tot_umi"
	gcor = tgs_cor(t(feat))
	diag(gcor) = NA

	dir_nm = scfigs_dir(gset_id, "gset_cors")
	if(!dir.exists(dir_nm)) {
		dir.create(dir_nm)
	}

	s_ids = unique(gset@gene_set)
	for(s_id in s_ids) {
		s_nms = names(gset@gene_set)[which(gset@gene_set==s_id)]
		s_nms = c(s_nms, "tot_umi")
		n = length(s_nms)
		fn = sprintf("%s/%s.png", dir_nm, s_id)
		png(fn, w=160+n*15, h=120+n*15)
		shades = colorRampPalette(c("darkblue", "blue","white", "red", "yellow"))(1000)
		breaks = seq(-0.5,0.5,l=1001)
		vs = pmin(pmax(gcor[s_nms,s_nms],-0.5),0.5)
		pheatmap::pheatmap(gcor[s_nms, s_nms], col=shades,breaks=breaks, cluster_rows= n > 2, cluster_cols= n > 2)
		dev.off()
	}
}
