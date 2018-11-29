#' plot a marker heat map give a metacell object
#'
#' @param mc_id id of metacell object ina scdb
#' @param gset_id markers to plot (e.g. using mcell_gset_from_mc_markers)
#' @param mat_id matrix object to us (default is the mc_id - assuming they use the same ID), if some cells in mc@mc are missing from the matrix, the function will generate an error
#' @param fig_fn file name for the figure (if null it will be call heat_marks in the fig directory)
#' @param lateral_gset_id markers to plot on top (e.g. cell cycle)
#' @param plot_cells by defulat this is TRUE and data is shown for single cells. If this is false, than metacells fold change values will be plotted
#'
#' @export

mcell_mc_plot_marks = function(mc_id, gset_id, mat_id = mc_id,
						fig_fn = NULL, lateral_gset_id = NULL,
						plot_cells = T, mc_ord = NULL)
{
	mcp_heatmap_height = get_param("mcp_heatmap_height")
	mcp_heatmap_width = get_param("mcp_heatmap_width")
	mcp_heatmap_ideal_umi = get_param("mcp_heatmap_ideal_umi")
	mcp_heatmap_fp_shades = colorRampPalette(get_param(ifelse(plot_cells, "mcp_heatmap_seq_shades", "mcp_heatmap_fp_shades")))(1000)
	mcp_heatmap_text_cex = get_param("mcp_heatmap_text_cex")
	mcp_heatmap_alt_side_text = get_param("mcp_heatmap_alt_side_text")
	mcp_heatmap_latgene_color = get_param("mcp_heatmap_latgene_color")

	mcp_heatmap_lwd = get_param("mcp_heatmap_lwd")
	mcp_heatmap_ord_cells_by_color_first = get_param("mcp_ord_cells_by_color_first")
	
	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
		stop("undefined meta cell object " , mc_id)
	}
	n_mc = ncol(mc@mc_fp)
	gset = scdb_gset(gset_id)
	if(is.null(gset)) {
		stop("undefined gset object when trying to plot markers, id " , gset_id)
	}
	scmat = scdb_mat(mat_id)
	if(is.null(scmat)) {
		stop("undefined mat object when trying to plot markers, id " , mat_id)
	}
	if(length(intersect(names(mc@mc), colnames(scmat@mat))) != length(names(mc@mc))) {
		stop("cells in meta cell are missing from provided matrix in mc_plot_marks")
	}
	if(is.null(mcp_heatmap_ideal_umi)) {
		mcp_heatmap_ideal_umi = quantile(colSums(scmat@mat), 0.25)
	}

	if(is.null(fig_fn)) {
		fig_fn = scfigs_fn(mc_id, ifelse(plot_cells, "cells_heat_marks", "mc_heat_marks"))
	}

	if (is.null(mc_ord)) {
		mc_ord = 1:ncol(mc@mc_fp)
	}
	cell_ord = names(mc@mc)[order(order(mc_ord)[mc@mc])]

	good_marks = intersect(names(gset@gene_set), rownames(mc@mc_fp))
	lateral_marks = c()
	if(!is.null(lateral_gset_id)) {
		lateral_gset = scdb_gset(lateral_gset_id)
		if(is.null(lateral_gset)) {
			stop("undefined lateral gset object when trying to plot markers, id " , lateral_gset_id)
		}
		lateral_marks = intersect(names(lateral_gset@gene_set), rownames(mc@mc_fp))
		good_marks = setdiff(good_marks, lateral_marks)
	}

	if(length(good_marks) < 2) {
		stop("Could not get >=2 markers to plot marker matrix")
	}

	gene_folds = mc@mc_fp

	if(!is.null(fig_fn)) {
	  if(is.null(mcp_heatmap_height)){
	    mcp_heatmap_height = 16*length(c(good_marks,lateral_marks)) + 250
	  }
	  if(is.null(mcp_heatmap_width)){
	    mcp_heatmap_width= max(min(3000,length(cell_ord)+200),800)
	  }
	  .plot_start(fig_fn, w=mcp_heatmap_width,h=mcp_heatmap_height);
	}
	layout(matrix(c(1,2),nrow=2),heights=c(mcp_heatmap_height, 100))

	top_marg=c(0,13,5,20)
	par(mar=top_marg)

	if(plot_cells) {
		mat = as.matrix(scmat@mat[c(good_marks, lateral_marks), cell_ord])
		totu = colSums(scmat@mat[, cell_ord])
		mat = t(t(mat)/totu)*mcp_heatmap_ideal_umi
		lus_1 = log2(1+7*mat)
		lus = apply(lus_1 - apply(lus_1, 1, median),2, function(x) pmax(x,0))
		if (length(cell_ord) < mcp_heatmap_width) {
			smooth_n = 1
		}
		smooth_n = max(2,ceiling(2*length(cell_ord)/mcp_heatmap_width))
		lus_smoo = t(apply(lus, 1, function(x) rollmean(x,smooth_n, fill=0)))

		
		image(t(lus_smoo), col=mcp_heatmap_fp_shades, xaxt='n', yaxt='n')

		vals_range = range(lus_smoo)
		
		cell_x = rep(NA, time=length(cell_ord))
		names(cell_x) = names(mc@mc)
		cell_x[cell_ord] = 1:length(cell_ord)
		cl_x = tapply(cell_x, mc@mc, mean)/length(cell_ord)
		mtext(1:n_mc, side = 3, at=cl_x, las=2, line = 2, cex=mcp_heatmap_text_cex)
		cl_x_b = tapply(cell_x, mc@mc, max)/length(cell_ord)
		abline(v=cl_x_b, lwd=mcp_heatmap_lwd)
	} else {
		mat = log2(gene_folds[c(good_marks, lateral_marks), mc_ord])
		mat = pmax(pmin(mat,3),-3)
		image(t(mat), col=mcp_heatmap_fp_shades, xaxt='n', yaxt='n', zlim=c(-3,3))
		mtext(mc_ord, side = 3, at=seq(0,1,l=n_mc), las=2, line = 2, cex=mcp_heatmap_text_cex)
		vals_range = c(-3, 3)
	}

	all_marks = c(good_marks, lateral_marks)
	g_n = length(all_marks)
	gene_cols = rep("black", g_n)
	if(length(lateral_marks) != 0) {
		gene_cols[(length(good_marks)+1):length(all_marks)] = mcp_heatmap_latgene_color
	}

	if(mcp_heatmap_alt_side_text) {
		odd = seq(1,g_n,2)
		even = seq(2,g_n,2)
		mtext(substr(all_marks[odd],1,8),
			at=seq(0,1,length.out=g_n)[odd],
			side=2, las=2, cex=mcp_heatmap_text_cex, col=gene_cols[odd])
		mtext(substr(all_marks[even],1,8),
			at=seq(0,1,length.out=g_n)[even],
			side=4, las=2, cex=mcp_heatmap_text_cex, col=gene_cols[even])
	} else {
		mtext(substr(all_marks,1,8),
			at=seq(0,1,length.out=g_n),
			side=2, las=2, cex=mcp_heatmap_text_cex, col=gene_cols)
		mtext(substr(all_marks,1,8),
			at=seq(0,1,length.out=g_n),
			side=4, las=2, cex=mcp_heatmap_text_cex, col=gene_cols)
	}

	lower_marg=c(5,13,0,20)
	par(mar=lower_marg)
	if(plot_cells) {
		image(as.matrix(1:length(cell_ord),nrow=1), col=mc@colors[mc@mc[cell_ord]], xaxt='n', yaxt='n')
		mtext(1:n_mc, side = 1, at=cl_x, las=2, line = 2, cex=mcp_heatmap_text_cex)
	} else {
		image(as.matrix(1:n_mc,nrow=1), col=mc@colors[mc_ord], yaxt='n', xaxt='n')
		mtext(mc_ord, side = 1, at=seq(0,1,l=n_mc), las=2, line = 2, cex=mcp_heatmap_text_cex)
	}
	dev.off()
	
	plot_color_bar(seq(vals_range[1], vals_range[2], length=1000), mcp_heatmap_fp_shades, gsub(".png", "_leg.png", fig_fn), show_vals_ind=c(1,500, 1000))
	
	write.table(all_marks, gsub(".png", "_gene_names.txt", fig_fn), quote=F, sep="\t")
	
}

#' PLot a series of heat maps to describe metacell groups
#'
#' @param mc_id metacell id in scdb
#' @param mat_id matrix object
#' @param lateral_gset_id genes to be marked as lateral
#'
#' @export


mcell_mc_plot_subheats = function(mc_id, mat_id, lateral_gset_id)
{
	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
		stop("undefined meta cell object " , mc_id)
	}
#cluster MCs into K groups by lfp
#for each group, create a temp subMC
#run the plotter
#delete the temp subMC
	tmp_mc_id = sprintf("%s_tmpsub", mc_id)
#	mc_hc = hclust(
#	for()
#		mcell_add_sub_mc(tmp_mc_id, mc_id, which(mc_clust==cl_i))
#		mcell_mc_plot_marks(tmp_mc_id, gset_id = tmp_gset_id,
#							mat_id = mat_id, fig_fn,
#							lateral_gset_id = lateral_gset_id)
#		scdb_del_mc(tmp_mc_id)
}
