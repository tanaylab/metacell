#' Plot mc+cell graph using pre-defined mc colorization
#'
#'

mcell_mc2d_plot = function(mc2d_id)
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
	fig_nm = scfigs_fn(mc2d_id, "2d_proj")
	png(fig_nm, width = mcp_2d_width, height = mcp_2d_height)
	cols = mc@colors
	cols[is.na(cols)] = "gray"
	plot(mc2d@sc_x, mc2d@sc_y, cex=mcp_2d_cex, pch=19, col=cols[mc@mc])
	fr = mc2d@graph$mc1
	to = mc2d@graph$mc2
	segments(mc2d@mc_x[fr], mc2d@mc_y[fr], mc2d@mc_x[to], mc2d@mc_y[to]) 
	points(mc2d@mc_x, mc2d@mc_y, cex= 3*mcp_2d_cex, col="black", pch=21, bg=cols)
	text(mc2d@mc_x, mc2d@mc_y, 1:length(mc2d@mc_x), cex=mcp_2d_cex)

	if(nrow(mc@color_key)!=0 & mcp_2d_plot_key) {
		key = mc@color_key
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
