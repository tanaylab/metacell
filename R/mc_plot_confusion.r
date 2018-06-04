#' plot a metacel confusion matrix
#'
#' @param mc_id id of metacell object ina scdb
#' @param graph_id cell to cell similarity graph
#' @param coc_id coclustering object to be used as the graph. If this is not null, graph_id must be null (and vice versa)
#' @param use_orig_order TRUE if you want to preserve the metacell order in the figu. Flase by default
#' @param mc_order defining an order for the MCs. If this is null and use_orig_order is false, the ordering will be based on hclust of the confusion matrix
#'
#' @export

mcell_mc_plot_confusion = function(mc_id, graph_id, coc_id = NULL,
							use_orig_order=F, mc_order =NULL, fig_fn = NULL)
{
	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
		stop("undefined meta cell object " , mc_id)
	}
	if(is.null(fig_fn)) {
		fig_fn = scfigs_fn(mc_id, sprintf("graph%s_confusion", graph_id))
	}

	if(!is.null(graph_id)) {
		if(!is.null(coc_id)) {
			stop("cannot specify both a graph and coclust graph when plotting confusion")
		}	
		cgraph = scdb_cgraph(graph_id)
		if(is.null(cgraph)) {
			stop("undefined cgraph object when trying to plot confusion, id " , graph_id)
		}
		max_deg = nrow(cgraph@edges)
		confu = mcell_mc_confusion_mat(mc_id, graph_id, max_deg, ignore_mismatch=T)
	} else if(!is.null(coc_id)) {
		coc = scdb_coclust(coc_id)
		if(is.null(coc)) {
			stop("undefined coclust object when trying to plot confusion, id " , coc_id)
		}
		max_deg = median(table(mc@mc))/2
		confu = mcell_mc_coclust_confusion_mat(mc_id, coc_id=coc_id, K=max_deg, ignore_mismatch=T, alpha=2)
	}
	r_confu = rowSums(confu)
	c_confu = colSums(confu)
	norm = r_confu %*% t(c_confu)
	confu_n = confu/norm

	colors = mc@colors

	if(!use_orig_order) {
		if(is.null(mc_order)) {
			hc = mcell_mc_hclust_confu(mc_id, graph_id=NULL, confu)
			mc_order = hc$order
		}
		confu = confu[mc_order, mc_order]
		confu_n = confu_n[mc_order, mc_order]
		colors = colors[mc_order]
		colnames(confu_n) = (1:ncol(confu_n))[mc_order]
		rownames(confu_n) = (1:ncol(confu_n))[mc_order]
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
}
