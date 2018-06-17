#' plot a "gel" like diagram showing expression of a gene of interest over metacells that are classified into types
#'
#' @param mc_id id of metacell object in scdb
#' @param gene_nm gene of interest
#' @param reorder_preset list of mc colors to define the "lanes"
#' @param reorder should metacell types be sorted by gene intesity?
#'
#' @export

mcell_mc_plot_vgels = function(mc_id, gene_nms,
							reorder_preset=NULL, reorder=F, fig_fn = NULL,
							lane_w=50, height=350)
{
	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
		stop("undefined meta cell object " , mc_id)
	}

	fdir = scfigs_dir(mc_id, "gels")
	for(gene_nm in gene_nms) {
		message("iter ", gene_nm)
		if(! gene_nm %in% rownames(mc@mc_fp)) {
			message("gene ", gene_nm, " is not covered/missing from mc footprint")
			next
		}
		fig_fn = sprintf("%s/%s.png", fdir, gene_nm)

		lfp = log2(mc@mc_fp)[gene_nm,]

		if(is.null(reorder_preset)) {
			if(reorder) {
				type_lfp = tapply(lfp, mc@colors, median)
				col_order = names(sort(type_lfp))
			} else {
				ord_lfp = tapply(1:length(lfp), mc@colors, median)
				col_order = names(sort(ord_lfp))
			}
		} else {
			col_order = reorder_preset
		}

		png(fig_fn, width=lane_w*length(col_order), height=height)
#		layout(matrix(c(1,2),nrow=2), heights=c(450,30))
		par(mar = c(0,5,4,3))

		max_lfp = max(lfp)
		min_lfp = min(min(lfp),0)
		max_lfp = max_lfp + (max_lfp-min_lfp)*0.02
		ybot = min_lfp - (max_lfp-min_lfp)/25
		ybot_up = min_lfp - (max_lfp-min_lfp)/70
		min_y = min_lfp - (max_lfp-min_lfp)/80

		plot(0,0, type='n', xlim=c(0,length(col_order)), ylim=c(min_y, max_lfp), xaxt='n', ylab=NA, xlab=NA)
		step = 0
		for(c in col_order) {
			shades = colorRampPalette(c("white", c))(402)[201:402]
			f_mc = mc@colors==c
			x = rank(runif(sum(f_mc))/1000+lfp[f_mc])/sum(f_mc)
			points(step+x, lfp[f_mc], pch=21, bg=shades[1+floor(x*200)], cex=1.5,lwd=0.5)
			step = step+1
		}
		grids = seq(0,length(col_order),l=length(col_order)+1)
		abline(v=grids, lwd=0.5)
		rect(xleft=grids[-1]-1, ybottom=ybot, xright=grids[-1],ytop=ybot_up, col=col_order)
#		image(as.matrix(1:length(col_order),nrow=1), col=col_order, yaxt='n', xaxt='n')
		dev.off()	
	}
}
