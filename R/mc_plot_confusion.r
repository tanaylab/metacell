#' plot a metacel confusion matrix
#'
#' @param mc_id id of metacell object ina scdb
#' @param graph_id cell to cell similarity graph
#' @param order_hc if this is true, the metacells will be reordered according to an hclust of the confusion matrix
#'
#' @export

mcell_mc_plot_confusion = function(mc_id, graph_id, order_hc=F, fig_fn = NULL)
{
	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
		stop("undefined meta cell object " , mc_id)
	}
	cgraph = scdb_cgraph(graph_id)
	if(is.null(cgraph)) {
		stop("undefined cgraph object when trying to plot confusion, id " , graph_id)
	}
	if(is.null(fig_fn)) {
		fig_fn = scfigs_fn(mc_id, sprintf("graph%s_confusion", graph_id))
	}
	
	max_deg = nrow(cgraph@edges)
	confu = mcell_mc_confusion_mat(mc_id, graph_id, max_deg, ignore_mismatch=T)
	r_confu = rowSums(confu)
	c_confu = colSums(confu)
	norm = r_confu %*% t(c_confu)
	confu_n = confu/norm

	colors = mc@colors
	if(order_hc) {
#		hc = hclust(dist(cor(log2(1+confu))),"ward.D2")
		confu_nodiag = confu_n
		diag(confu_nodiag) = 0
		confu_n = pmin(confu_n, max(confu_nodiag))
		confu_n = pmin(confu_n, quantile(confu_n, 1-3/nrow(confu_n)))
		epsilon = quantile(confu_n[confu_n!=0],0.02)
		hc = hclust(as.dist(-log10(epsilon+confu_n)),"average")
		confu = confu[hc$order, hc$order]
		confu_n = confu_n[hc$order, hc$order]
		colors = colors[hc$order]
	}
	colors[is.na(colors)] = "gray"

	shades = colorRampPalette(c("white", "pink", "red", "black", "brown", "orange"))
	png(fig_fn, w=800, h=800)
	layout(matrix(c(1,4,2,3),nrow=2),heights=c(800, 50), width=c(50,800))

	tl_marg=c(0,2,5,0)
	par(mar=tl_marg)
	n_mc = ncol(mc@mc_fp)
	image(t(as.matrix(1:n_mc,nrow=1)), col=colors, yaxt='n', xaxt='n')

	top_marg=c(0,0,5,5)
	par(mar=top_marg)
	log_scale = F
	if(log_scale) {
		image(log2(1+confu),col=shades(1000),xaxt='n', yaxt='n')
	} else {
		confu_nodiag = confu_n
		diag(confu_nodiag) = 0
		confu_n = pmin(confu_n, max(confu_nodiag))
		confu_n = pmin(confu_n, quantile(confu_n, 1-3/nrow(confu_n)))
		image(confu_n,col=shades(1000),xaxt='n', yaxt='n')
	}

	lower_marg=c(3,0,0,5)
	par(mar=lower_marg)
	image(as.matrix(1:n_mc,nrow=1), col=colors, yaxt='n', xaxt='n')
	dev.off()

	tree_fig_fn = sub("png", "hc.png", fig_fn)
	png(tree_fig_fn, w = 1000, h=400)
	plot(hc, cex=0.1)
	grid()
	dev.off()

	return(hc)
}
