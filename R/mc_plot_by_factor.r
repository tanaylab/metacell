#' plot a heatmap of number of cells per metacell and metadata factor (e.g. patient, condition, sample etc.)
#'
#' @param mc_id id of metacell object ina scdb
#' @param mat_id matrix object to us (default is the mc_id - assuming they use the same ID), if some cells in the mc object are missing from the matrix, the function will generate an error
#' @param fig_fn file name for the figure. Color bars are generated in separate files with names starting with fig_fn. If null figures will be created under <mc_id>.mc_comp_by_<meta_field> directory, with name composed of the plotting parameters.
#' @param meta_field metdata column name used to break down cells on metacells
#' @param meta_field_annotate_by metadata column names describing the meta_field (should be a one-to-one match between each meta_field_value and the values in these columns). Expecting a yaml parameter named mcp_metadata_annot_colors which is a list mapping meta_field_annotate_by values to list of value-color mapping.
#' @param meta_field_min_count Minimal number of cells per meta_field value 
#' @param norm_by_factor Normalize cell counts by factor (default) or by metacell
#' @param hclust_mcs If true, re-order metacells by clustering the heatmap. Otherwise, use original metacells ordering (default)
#' @param filter_values vector of meta_field values to display
#'

mcell_mc_plot_by_factor = function(mc_id, meta_field, mat_id = mc_id, fig_fn = NULL, meta_field_annotate_by=NULL, meta_field_min_count=0, norm_by_factor=T, hclust_mcs=F, ord_mcs = NULL, filter_values=NULL)
{
	mcp_heatmap_width = get_param("mcp_heatmap_width")
	mcp_heatmap_shades = get_param("mcp_heatmap_seq_shades")
   mcp_metadata_annot_colors = get_param("mcp_metadata_annot_colors")
   mcp_heatmap_text_cex = get_param("mcp_heatmap_text_cex")

	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
	  stop(sprintf("undefined metacell object when trying to plot metacells by %s, id = %s", meta_field , mc_id))
	}

	scmat = scdb_mat(mat_id)
	if(is.null(scmat)) {
		stop(sprintf("undefined mat object when trying to plot metacells by %s, id = %s", meta_field, mat_id))
	}
	if(length(intersect(names(mc@mc), colnames(scmat@mat))) != length(names(mc@mc))) {
		stop("cells in meta cell are missing from provided matrix in mcell_mc_plot_by_factor")
	}

	if (!(meta_field %in% colnames(scmat@cell_metadata))) {
		stop(sprintf("Column named %s not found in metadata of mat_id %s", meta_field, mat_id))
	}
	mc_t = table(mc@mc, scmat@cell_metadata[names(mc@mc), meta_field])
	mc_t = mc_t[, colSums(mc_t) >= meta_field_min_count]
	if (!is.null(filter_values)) {
		mc_t = mc_t[, colnames(mc_t) %in% filter_values]
	}
	if (ncol(mc_t) == 0) {
		stop(sprintf("nothing to plot in mcell_mc_plot_by_factor by %s when requiring at least %d cells %s for factor %s", meta_field, meta_field_min_count, ifelse(is.null(filter_values), "", paste0(" and these values: ", paste0(filter_values, collapse=", ")))))
	}
	
	if (norm_by_factor) {
		mc_t_n = t(t(mc_t) / colSums(mc_t))
	}
	else {
		mc_t_n = mc_t / rowSums(mc_t)
	}

	if(is.null(fig_fn)) {
		fig_fn = scfigs_fn(mc_id, 
											 sprintf("mc_comp_by_%d_%s_norm%s%s", ncol(mc_t), meta_field, ifelse(norm_by_factor, meta_field, "MC"), ifelse(hclust_mcs, "_hclust_mcs", "")),
											 sprintf("%s/%s.mc_comp_by_%s", .scfigs_base, mc_id, meta_field))
	}
	
	.plot_start(fig_fn, w=mcp_heatmap_width,h=100 + ncol(mc_t) * 15)

	n_ann = 1 + ifelse(is.null(meta_field_annotate_by), 0, length(meta_field_annotate_by))
  layout(matrix(c(rep(n_ann+3, n_ann), n_ann+1, 1:n_ann, n_ann+2), nrow=2, byrow=T), heights = c(2, ncol(mc_t)), widths = c(rep(1, n_ann), n_ann * 16))

	par(mar=c(1,1,1,10))

	chc = hclust(dist(t(mc_t_n)), method="ward.D2")
	if(!is.null(ord_mcs)) {
		rord = ord_mcs
	} else if(hclust_mcs) {
		rhc = hclust(dist(mc_t_n), method="ward.D2")
		rord = rhc$order
	} else {
		rord = 1:nrow(mc_t_n)
	}

	mf_ann = data.frame(row.names=colnames(mc_t), n_cells=colSums(mc_t))
	if (!is.null(meta_field_annotate_by)) {
		missing_fields = setdiff(meta_field_annotate_by, colnames(scmat@cell_metadata))
		if (length(missing_fields) > 0) {
			stop(sprintf("Columns %s not found in metadata of mat_id %s", paste0(missing_fields, collapse=", "), mat_id))
		}
	  uniq_md = unique(scmat@cell_metadata[names(mc@mc), c(meta_field, meta_field_annotate_by)])
	  if (nrow(uniq_md) != length(unique(uniq_md[, meta_field]))) {
	    stop(sprintf("annotation fields (%s) are not unique for metadata field %s", paste0(meta_field_annotate_by, collapse=", "), meta_field))
	  }
	  rownames(uniq_md) = uniq_md[, meta_field]
	  mf_ann = cbind(mf_ann, uniq_md[rownames(mf_ann), meta_field_annotate_by])
	  mf_ann[is.na(mf_ann)] = "N/A"
	  mf_ann[mf_ann == ""] = "N/A"

	  mf_ann = mf_ann[chc$order, ]
	}

	par(mar=c(10,0.5,0.5,0.5))

	for (annot in meta_field_annotate_by) {
	  annot_vals = unlist(mcp_metadata_annot_colors[[annot]])
	  image(t(1:nrow(mf_ann)), col=annot_vals[mf_ann[, annot]], xaxt='n', yaxt='n', xlab="", ylab="")
	  mtext(annot, 1, line=0.5, cex = mcp_heatmap_text_cex, las=2)
	}

	n_cells_zlim = c(0, 100 * ceiling(max(mf_ann$n_cells)/100))
	par(mar=c(10,0.5,0.5,0.5))
	image(t(mf_ann$n_cells), zlim=n_cells_zlim, col=colorRampPalette(RColorBrewer::brewer.pal(n=9, 'Blues'))(101), xaxt='n', yaxt='n', xlab="", ylab="")
	mtext("#cells", 1, line=0.5, cex = mcp_heatmap_text_cex, las=2)

	par(mar=c(0.5, 0.5, 0.5, 20))
	image(as.matrix(1:nrow(mc_t)), col=mc@colors[rord], xaxt='n', yaxt='n', xlab="", ylab="")

	par(mar=c(10, 0.5, 0.5, 20))
	mc_t_n = mc_t_n[rord, chc$order]
	#mc_t_enr = mc_t_enr[, chc$order]

	max_p = ceiling(100 * quantile(mc_t_n, 0.99))/100
	image(pmin(mc_t_n, max_p), col=colorRampPalette(mcp_heatmap_shades)(101), zlim=c(0, max_p), xaxt='n', yaxt='n', xlab="", ylab="")
	mtext(text=rownames(mc_t_n), at=seq(0, 1, length=nrow(mc_t_n)), side=1, line=0.5, cex = mcp_heatmap_text_cex/2, las=2)
	mtext(colnames(mc_t_n), at=seq(0, 1, length=ncol(mc_t_n)), side=4, line=0.5, cex = mcp_heatmap_text_cex, las=2)
	abline(h=seq(0.5/ncol(mc_t_n), 1, by=1/(ncol(mc_t_n)-1)), col="grey80")

  dev.off()

	# plot color bars
	plot_color_bar(vals=seq(0, max_p, length=101), cols=colorRampPalette(mcp_heatmap_shades)(101), title="fraction", show_vals_ind=c(1, 101), fig_fn=gsub(".png", "_color_bar_main.png", fig_fn))

	plot_color_bar(vals=seq(n_cells_zlim[1], n_cells_zlim[2], length=101), cols=colorRampPalette(RColorBrewer::brewer.pal(n=9, 'Blues'))(101), title="# cells", show_vals_ind=c(1, 101), fig_fn=gsub(".png", "_color_bar_n_cells.png", fig_fn))

	for (annot in meta_field_annotate_by) {
		if (is.null(mcp_metadata_annot_colors[[annot]])) {
			stop(sprintf("Missing value-color mapping entry for %s in parameter mcp_metadata_annot_colors", annot))
		}
	  annot_vals = rev(unlist(mcp_metadata_annot_colors[[annot]]))
	  annot_vals = annot_vals[ names(annot_vals) %in% mf_ann[, annot] ]
	  plot_color_bar(vals=names(annot_vals), cols=annot_vals, title=annot, fig_fn=gsub(".png", paste0("_color_bar_", annot, ".png"), fig_fn))
	}
	
}
