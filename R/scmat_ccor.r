#' Analyze cell cell cor
#'
#' 
#'
#' @param mat_id matrix to use
#' @param gset_id features to use for computing cors
#' @param n_dsamp number of umis to for downsampling. If this is -1, no dsamp will be performed
#'
#' @export
#'

mcell_scmat_plot_cmp_kcor = function(mat_id, gset_id)
{
	feat = gset_get_feat_mat(gset_id, mat_id, downsamp=T)
	k_nonz_exp = get_param("scm_k_nonz_exp")
#force given dsamp?
	feat = log2(1+k_nonz_exp*as.matrix(feat))

	raw_cor = tgs_cor(feat)
	
#distribution of k'th top cor	
	fn = scfigs_fn(mat_id, sprintf("k_cor_cmp_%s", gset_id))
	png(fn, w=1200, h=600)
	layout(matrix(c(1,5,2,6,3,7,4,8), ncol=4))
	k_titer = apply(raw_cor, 1, max)
	diag(raw_cor) = 0
	min_cor = min(raw_cor) - 0.03
	max_cor = max(raw_cor) + 0.03
	max_lk = min(floor(log2(ncol(feat))), 11)
	ref_cor = apply(raw_cor, 1, function(x) { -sort(-x, partial=4)[4] })
	for(k in 2**(4:max_lk)) {
		k_cor = apply(raw_cor, 1, function(x) { -sort(-x, partial=k)[k] })
		plot(ref_cor, k_cor, pch=19, cex=1, 
					main=sprintf("comparing K=%d",k), 
					xlim=c(min_cor, max_cor), ylim=c(min_cor, max_cor), 
					xlab="knn k=4", ylab=sprintf("knn k= %d", k))
		grid()
		abline(b=1,a=0, lty=2)
	}
	dev.off()
	rank_out = apply(raw_cor, 1, rank)
	rank_in = t(rank_out)
	f = rank_out<200
	fn = scfigs_fn(mat_id, "cor_in_out") 
	png(fn, w=800, h=800)
	cor_shades = colorRampPalette(c("darkblue", "blue", "black", "red", "yellow"))(201)
	smoothScatter(as.vector(rank_out)[f], log2(1+as.vector(rank_in)[f]), 
				nrpoints = ncol(rank_out)*2, cex=0.4, pch=19, 
				colramp = colorRampPalette(c("white", "blue", "orange")))
#				col = cor_shades[ceiling(100*(1+as.vector(raw_cor)[f]))])
	dev.off()
	return(list(rank_out, rank_in))
}
