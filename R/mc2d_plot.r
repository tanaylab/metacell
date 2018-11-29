#' Plot mc+cell graph using pre-defined mc colorization
#'
#' @param mc2d_id mc2d object to plot
#' @param legend_pos position of legend
#' @param plot_edges plot edges between metacells (true by default)
#'
#' @export
mcell_mc2d_plot = function(mc2d_id, legend_pos="topleft", plot_edges=T)
{
	mcp_2d_height = get_param("mcell_mc2d_height")
	mcp_2d_width = get_param("mcell_mc2d_width")
	mcp_2d_plot_key = get_param("mcell_mc2d_plot_key")
	mcp_2d_cex = get_param("mcell_mc2d_cex")
	mcp_2d_legend_cex = get_param("mcell_mc2d_legend_cex")

	mc2d = scdb_mc2d(mc2d_id)
	if(is.null(mc2d)) {
		stop("missing mc2d when trying to plot, id ", mc2d_id)
	}
	mc = scdb_mc(mc2d@mc_id)
	if(is.null(mc)) {
		stop("missing mc in mc2d object, id was, ", mc2d@mc_id)
	}
	fig_nm = scfigs_fn(mc2d_id, ifelse(plot_edges, "2d_graph_proj", "2d_proj"))
	.plot_start(fig_nm, w=mcp_2d_width, h = mcp_2d_height)
	#png(fig_nm, width = mcp_2d_width, height = mcp_2d_height)
	cols = mc@colors
	cols[is.na(cols)] = "gray"
	plot(mc2d@sc_x, mc2d@sc_y, pch=19, col=cols[mc@mc[names(mc2d@sc_x)]])
	fr = mc2d@graph$mc1
	to = mc2d@graph$mc2
	if (plot_edges) {
		segments(mc2d@mc_x[fr], mc2d@mc_y[fr], mc2d@mc_x[to], mc2d@mc_y[to])
	}
	points(mc2d@mc_x, mc2d@mc_y, cex= 3*mcp_2d_cex, col="black", pch=21, bg=cols)
	text(mc2d@mc_x, mc2d@mc_y, 1:length(mc2d@mc_x), cex=mcp_2d_cex)

	if(nrow(mc@color_key)!=0 & mcp_2d_plot_key) {
		key = mc@color_key[ mc@color_key$color %in% mc@colors, ]
#group	gene	color	priority	T_fold
		gmark = tapply(key$gene, key$group, paste, collapse=", ")
		gcol = unique(data.frame(col=key$color, group=key$group))
		rownames(gcol) = gcol$group
		gmark = gmark[order(names(gmark))]
		legend(legend_pos,
				legend=gsub("_", " ", paste0(names(gmark), ": ", gmark)),
				pch=19, cex=mcp_2d_legend_cex,
				col=as.character(gcol[names(gmark), 'col']), bty='n')
	}

	dev.off()
}


#' Plot mc+cells using pre-defined mc colorization, breakdown by given metadata field (e.g. patient)
#'
#' @param mc2d_id mc2d object to use for plot
#' @param mat_id mat object matching mc2d_id that contains the cells metadata information
#' @param meta_field field name (in mat cell_metadata slot) to split cells by
#' @param single_plot output all panels in a single plot or plot per panel (T)
#' @param filter_values to filter meta_field values by (NULL)
#' @param filter_name name to add to plots (NULL)
#' @param ncols number of panels in column (if single_plot is true), automatically calculate if NULL
#' @param neto_points plot without a box and a title (relevant if single_plot=F)
#'
#' @export
mcell_mc2d_plot_by_factor = function(mc2d_id, mat_id, meta_field, single_plot = T, filter_values = NULL, filter_name = NULL, ncols=NULL, neto_points=F)
{
  mcp_2d_height = get_param("mcell_mc2d_height")
  mcp_2d_width = get_param("mcell_mc2d_width")
  mcp_2d_plot_key = get_param("mcell_mc2d_plot_key")
  mcp_2d_cex = 2 * get_param("mcell_mc2d_cex")
  mcp_2d_legend_cex = get_param("mcell_mc2d_legend_cex")

  bg_col="grey90"

  mc2d = scdb_mc2d(mc2d_id)
  if(is.null(mc2d)) {
    stop("missing mc2d when trying to plot, id ", mc2d_id)
  }
  mc = scdb_mc(mc2d@mc_id)
  if(is.null(mc)) {
    stop("missing mc in mc2d object, id was, ", mc2d@mc_id)
  }

  mat = scdb_mat(mat_id)
  if (is.null(mat)) {
    stop(sprintf("missing mat (id = %s) for metadata info when plotting %s", mat_id, mc2d_id))
  }

  if (any(mat@cells != mc@cell_names)) {
    stop(sprintf("cells mismatch between mc2d mc (id = %s) and mat (id = %s) objects", mc2d_id, mat_id))
  }

  c_by_f = split(names(mc@mc), mat@cell_metadata[names(mc@mc), meta_field])
  if (is.null(filter_values)) {
    filter_values = names(c_by_f)
  }
  else {
    c_by_f = c_by_f[names(c_by_f) %in% filter_values]
  }

  cols = mc@colors
  cols[is.na(cols)] = "white"

  if (single_plot) {
    ny = ifelse(is.null(ncols), floor(sqrt(length(c_by_f))), ncols)
    nx = ceiling((length(c_by_f)/ny))

    fig_nm = scfigs_fn(mc2d_id, sprintf("2d_proj_%sall", ifelse(is.null(filter_name), "", paste0(filter_name, "_"))), sprintf("%s/%s.by_%s", .scfigs_base, mc2d_id, meta_field))
    .plot_start(fig_nm, w=mcp_2d_width, h=mcp_2d_width / ny * nx)

    layout(matrix(1:(nx*ny), nx, ny, byrow=T))
    par(mar=c(0.5,0.5,3,0.5))
  }

  for (meta_field_v in filter_values) {
    ccells = c_by_f[[meta_field_v]]

    if (!single_plot) {
      fig_nm = scfigs_fn(mc2d_id, sprintf("2d_proj_%s", meta_field_v, ifelse(is.null(filter_name), "", paste0(filter_name, "_"))), sprintf("%s/%s.by_%s", .scfigs_base, mc2d_id, meta_field))
      .plot_start(fig_nm, w=mcp_2d_width, h=mcp_2d_height)
      par(mar=c(0.5, 0.5, ifelse(neto_points, 0.5, 3), 0.5))
    }

    #col=cols[mc@mc]
    plot(mc2d@sc_x, mc2d@sc_y, cex=mcp_2d_cex, pch=21, col=bg_col, bg=bg_col, xlab="", ylab="", xaxt='n', yaxt='n', bty=ifelse(!single_plot & neto_points, 'n', 'o'))
    points(mc2d@sc_x[ccells], mc2d@sc_y[ccells], cex= mcp_2d_cex, lwd=0.5, pch=19, col=cols[mc@mc[ccells]])
    #text(mc2d@mc_x, mc2d@mc_y, 1:length(mc2d@mc_x), cex=mcp_2d_cex)

    if (single_plot || !neto_points) {
    	title(main=sprintf("%s (%d)", meta_field_v, length(ccells)), cex.main=1.2*mcp_2d_legend_cex)
    }
    if (!single_plot) {

      if(nrow(mc@color_key) != 0 & mcp_2d_plot_key & !neto_points) {
        key = mc@color_key[ mc@color_key$color %in% mc@colors, ]
        #group	gene	color	priority	T_fold
        gmark = tapply(key$gene, key$group, paste, collapse=", ")
        gcol = unique(data.frame(col=key$color, group=key$group))
        rownames(gcol) = gcol$group
        gmark = gmark[order(names(gmark))]
        legend("topleft",
               legend=gsub("_", " ", paste0(names(gmark), ": ", gmark)),
               pch=19, cex=mcp_2d_legend_cex,
               col=as.character(gcol[names(gmark), 'col']), bty='n')
      }

      dev.off()
    }
  }

  if (single_plot) {
    dev.off()
  }

}

#' Plot the (log2) metacell footprint value of the selected gene on the 2d projection
#'
#' @param mc2d_id mc2d object to use for plot
#' @param gene gene name to plot
#' @param show_mc_ids plot metacell ids (false by default)
#' @param show_legend plot color bar legend (true by default)
#' @param neto_points do not plot box, title and legend
#'
#' @export
#'
mcell_mc2d_plot_gene = function(mc2d_id, gene, show_mc_ids=F, show_legend=T, neto_points=F)
{
	height = get_param("mcell_mc2d_gene_height")
	width = get_param("mcell_mc2d_gene_width")
	mc_cex = get_param("mcell_mc2d_gene_mc_cex")
	sc_cex = get_param("mcell_mc2d_gene_cell_cex")
	colspec = get_param("mcell_mc2d_gene_shades")
	max_lfp = get_param("mcell_mc2d_gene_max_lfp")
	
	mc2d = scdb_mc2d(mc2d_id)
	if(is.null(mc2d)) {
		stop("missing mc2d when trying to plot, id ", mc2d_id)
	}
	mc = scdb_mc(mc2d@mc_id)
	if(is.null(mc)) {
		stop("missing mc in mc2d object, id was, ", mc2d@mc_id)
	}
	
	if (!(gene %in% rownames(mc@mc_fp))) {
		stop(sprintf("gene %s not found in mc object id %s mc_fp table", gene, mc2d@mc_id))
	}
	
	x = pmin(pmax(log2(mc@mc_fp[gene, ]), -max_lfp), max_lfp) + max_lfp
	shades = colorRampPalette(colspec)(200 * max_lfp + 1)
	mc_cols = shades[round(100 * x) + 1] 
	
	fig_nm = scfigs_fn(mc2d_id, gene, sprintf("%s/%s.genes", .scfigs_base, mc2d_id))
	.plot_start(fig_nm, w = width * ifelse(show_legend & !neto_points, 1.25, 1), h = height)
	if (show_legend & !neto_points) {
		layout(matrix(c(1,1:3), nrow=2, ncol=2), widths = c(4,1))
	}
	
	if (neto_points) {
		par(mar=c(1,1,1,1))
	} else {
		par(mar=c(4,4,4,1))
	}

	plot(mc2d@sc_x, mc2d@sc_y, pch=19, cex=sc_cex, col='grey', xlab="", ylab="", main=ifelse(neto_points, "", gene), cex.main=mc_cex, bty=ifelse(neto_points, 'n', 'o'), xaxt=ifelse(neto_points, 'n', 's'), yaxt=ifelse(neto_points, 'n', 's'))
	points(mc2d@mc_x, mc2d@mc_y, pch=21, bg=mc_cols, cex=mc_cex)
	
	if (show_mc_ids) {
		text(mc2d@mc_x, mc2d@mc_y, seq_along(mc2d@mc_y), cex=mc_cex * 0.5)
	}
	
	if (show_legend & !neto_points) {
		par(mar=c(4,1,4,1))
		plot_color_bar(seq(-max_lfp, max_lfp, l=length(shades)), shades, show_vals_ind=c(1, 100 * max_lfp + 1, 200 * max_lfp + 1))
	}
	
	dev.off()
}
