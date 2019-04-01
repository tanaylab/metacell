#' plot gene/feature statistics
#'
#' @param gstat_id id of gene state object
#' @param gset_id optional gene set, will be used to mark select genes if specified
#' @param fig_dir optional base dir (default is to use the std fig dir)
#'
#' @export
#'
mcell_plot_gstats = function(gstat_id, gset_id = NULL, fig_dir = NULL, max_vm=4)
{
	gstat = scdb_gstat(gstat_id)
	if(is.null(gstat)) {
		stop("MC-ERR non existing gstats in mcell_plot_gstats, id ", gstat_id)
	}
	if(!is.null(gset_id)) {
		gset = scdb_gset(gset_id)
		if(is.null(gset)) {
			stop("MC-ERR non existing gset id ", gset_id, " when calling plot gstat")
		}
		marks = names(gset@gene_set)

		set_cols = rep(RColorBrewer::brewer.pal(n=max(3, min(length(gset@set_names), 9)) , 'Set1'), times=length(gset@set_names))[1:length(gset@set_names)]
		names(set_cols) = gset@set_names

		cols = unlist(set_cols[gset@gene_set])

	}
	if(is.null(fig_dir)) {
		vm_fn = scfigs_fn(gstat_id, "varmin")
		szcor_fn = scfigs_fn(gstat_id, "szcor")
		top3_fn = scfigs_fn(gstat_id, "top3")
	} else {
		vm_fn = sprintf("%s/%s_varmin.png", fig_dir, gstat_id)
		szcor_fn = sprintf("%s/%s_szcor.png", fig_dir, gstat_id)
		top3_fn = sprintf("%s/%s_top3.png", fig_dir, gstat_id)
	}
	png(vm_fn, w=1200, h=1200, pointsize=7, res=300)
	vm = pmin(gstat$ds_log_varmean, ifelse(is.null(max_vm), max(gstat$ds_log_varmean), max_vm))
	names(vm) = rownames(gstat)
   plot(log2(gstat$ds_mean), gstat$ds_log_varmean, cex=0.8, pch=19,
		 ylim = c(min(vm), max(vm)+0.25),
       xlab = "log(downsampled mean)", ylab="log2(var/mean downsampled)");
	if(!is.null(gset_id)) {
		points(log2(gstat[marks,"ds_mean"]),
				 vm[marks],
				 cex=0.8, pch=19, col=cols);
	  legend("topleft", legend=names(set_cols), pch=19, cex=0.8, col=set_cols, bty='n')
	}
	dev.off()

	png(szcor_fn, w=1200, h=1200, pointsize=7, res=300)
   plot(log2(gstat$ds_mean), gstat$sz_cor, cex=0.8, pch=19,
       xlab = "log(downsampled mean)", ylab="sz correlation");
	if(!is.null(gset_id)) {
		points(log2(gstat[marks,"ds_mean"]),
				 gstat[marks, "sz_cor"],
				 cex=0.8, pch=19, col=cols);
	  legend("topleft", legend=names(set_cols), pch=19, cex=0.8, col=set_cols, bty='n')
	}
	dev.off()

	png(top3_fn, w=1200, h=1200, pointsize=7, res=300)
   plot(log2(gstat$ds_mean), log2(1+gstat$ds_top3), cex=0.8, pch=19,
       xlab = "log(downsampled mean)", ylab="log third highest umi (downsamp)");
	if(!is.null(gset_id)) {
		points(log2(gstat[marks,"ds_mean"]),
				 log2(1+gstat[marks, "ds_top3"]),
				 cex=0.8, pch=19, col=cols);
	  legend("topleft", legend=names(set_cols), pch=19, cex=0.8, col=set_cols, bty='n')
	}
	dev.off()
}
#
