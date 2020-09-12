#' plot a marker heat map give a metacell object
#'
#' @param mc_id id of metacell object ina scdb
#' @param gset_id markers to plot (e.g. using mcell_gset_from_mc_markers)
#' @param mat_id matrix object to us (default is the mc_id - assuming they use the same ID), if some cells in mc@mc are missing from the matrix, the function will generate an error
#' @param fig_fn file name for the figure (if null it will be call heat_marks in the fig directory)
#' @param lateral_gset_id markers to plot on top (e.g. cell cycle)
#' @param plot_cells by defulat this is TRUE and data is shown for single cells. If this is false, than metacells fold change values will be plotted
#' @param zero_median if this is true then for each gene, all values below the median will be plotted as "0". It can save some noise but will  make genes that are anit-enriched in some MC non visible.
#' @param add_genes A list of marker genes to be added manually
#' @param add_metadata name of metadata field to add below hte heat map
#' @param ext_metadata values of metadata per cell to add (a anamed vector, with some or all of the cell names specified)
#' @param md_level_colors what colors to use for the metadata levels, if present
#' @param focus_mcs option list of MC IDs in case only some of the model should be plotted
#' @param gene_list a list of gene names to replcace the markers from a gset_id
#'
#' @export

mcell_mc_plot_marks = function(mc_id, gset_id, mat_id = mc_id,
						fig_fn = NULL, lateral_gset_id = NULL,
						mc_ord = NULL,
						plot_cells = T,
						zero_median = T,
						add_genes = NULL,
						add_metadata=NULL,
						ext_metadata=NULL,
						md_level_colors = NULL,
						focus_mcs = NULL,
						gene_list = NULL,
						reorder_marks=F, fold_burn=3)
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
	if(!is.null(focus_mcs)) {
		mc = mc_restrict(mc, focus_mcs)
	}
	
	n_mc = ncol(mc@mc_fp)
	if(is.null(gene_list)) {
	gset = scdb_gset(gset_id)
		if(is.null(gset)) {
			stop("undefined gset object when trying to plot markers, id " , gset_id)
		}
		gene_list = names(gset@gene_set)
	}
	scmat = scdb_mat(mat_id)
	if(is.null(scmat)) {
		stop("undefined mat object when trying to plot markers, id " , mat_id)
	}
	if(length(intersect(names(mc@mc), colnames(scmat@mat))) != length(names(mc@mc))) {
		stop("cells in meta cell are missing from provided matrix in mc_plot_marks")
	}
	n_md_levels = 0
	plot_meta = F
	if(!is.null(add_metadata)) {
		if(!add_metadata %in% colnames(scmat@cell_metadata)) {
			stop("unknown meta data field " , add_metadata, " in plot marks")
		}
		n_md_levels = length(levels(as.factor(as.character(scmat@cell_metadata[names(mc@mc),add_metadata]))))
		message("will plot metadata with ", n_md_levels, " levels")
		if(n_md_levels > 50) {
			stop("meta data field " , add_metadata, " has ",n_md_levels, " levels, while maximum is 50 - try to bin it maybe?")
		}
		plot_meta = T
	}
	if(!is.null(ext_metadata)) {
		if(is.null(names(ext_metadata)) |
		  length(intersect(names(ext_metadata), names(mc@mc))) == 0) {
			stop("External metadata should be a list with names being (a subset) of the cells in the metacell object")
		}
		n_md_levels = length(levels(as.factor(as.character(ext_metadata))))
		plot_meta = T
	}
	if(is.null(mcp_heatmap_ideal_umi)) {
		mcp_heatmap_ideal_umi = quantile(colSums(scmat@mat), 0.25)
	}
	if(n_md_levels == 1) {
		message("only one metadata level, will supress plotting")
		plot_meta=F
	}

	if(is.null(fig_fn)) {
		fig_fn = scfigs_fn(mc_id, ifelse(plot_cells, "cells_heat_marks", "mc_heat_marks"))
	}

	good_marks = intersect(gene_list, rownames(mc@mc_fp))
	if(!is.null(focus_mcs)) {
		mc_ord = 1:length(focus_mcs)
		names(mc_ord) = focus_mcs
		cell_ord = names(mc@mc)[order(mc_ord[as.character(mc@mc)]+runif(length(mc@mc),max=0.02))]
		mc_ord = focus_mcs
	} else if(is.null(mc_ord)) {
#	   mc_ord = 1:ncol(mc@mc_fp)
		mc_ord = 1:length(unique(mc@mc))
		cell_ord = names(mc@mc)[order(order(mc_ord)[mc@mc]+runif(length(mc@mc),max=0.02))]
	} else {
		fps = mc@mc_fp[good_marks,mc_ord]
		if(reorder_marks) {
			good_marks = good_marks[order(apply(fps,1, which.max))]
		}
		cell_ord = names(mc@mc)[order(order(mc_ord)[mc@mc]+runif(length(mc@mc),max=0.02))]
	}

	lateral_marks = c()
	if(!is.null(lateral_gset_id)) {
		lateral_gset = scdb_gset(lateral_gset_id)
		if(is.null(lateral_gset)) {
			stop("undefined lateral gset object when trying to plot markers, id " , lateral_gset_id)
		}
		lateral_marks = intersect(names(lateral_gset@gene_set), rownames(mc@mc_fp))
		good_marks = setdiff(good_marks, lateral_marks)
	}
	if(!is.null(add_genes)) {
		good_marks = c(good_marks, intersect(add_genes, rownames(mc@mc_fp)))
	}

	if(length(good_marks) < 2) {
		stop("Could not get >=2 markers to plot marker matrix")
	}

	gene_folds = mc@mc_fp

	if(!is.null(fig_fn)) {
	  if(is.null(mcp_heatmap_height)){
	    mcp_heatmap_height = 16*length(c(good_marks,lateral_marks)) + 650 + n_md_levels*25
	  }
	  if(is.null(mcp_heatmap_width)){
	    mcp_heatmap_width= max(min(3000,length(cell_ord)+200),800)
	  }
	  message("setting fig h to ", mcp_heatmap_height, " md levels ", n_md_levels, " num of marks " , length(good_marks))
	  png(fig_fn, w=mcp_heatmap_width,h=mcp_heatmap_height);
	}
	if(plot_meta) {
		layout(matrix(c(1,2,3),nrow=3),heights=c(mcp_heatmap_height, 50, 50+n_md_levels*25))
	} else {
		layout(matrix(c(1,2),nrow=2),heights=c(mcp_heatmap_height, 100))
	}

	top_marg=c(0,13,5,20)
	par(mar=top_marg)
	if(plot_cells) {
		mat = as.matrix(scmat@mat[c(good_marks, lateral_marks),names(mc@mc)])
		mat = mat[,cell_ord]
		totu = colSums(scmat@mat[,cell_ord])
		mat = t(t(mat)/totu)*mcp_heatmap_ideal_umi
		lus_1 = log2(1+7*mat)
		if(zero_median) {
			lus = apply(lus_1 - apply(lus_1, 1, median),2, function(x) pmax(x,0))
		} else {
			lus = lus_1 - apply(lus_1, 1, median)
			lus = pmin(pmax(lus,-3),3)
		}
		if (length(cell_ord) < mcp_heatmap_width) {
			smooth_n = 1
		}
		smooth_n = max(2,ceiling(2*length(cell_ord)/mcp_heatmap_width))
		lus_smoo = t(apply(lus, 1, function(x) rollmean(x,smooth_n, fill=0)))
		if(reorder_marks) {
			hc_g = hclust(tgs_dist(lus_smoo[good_marks,]),"ward.D2")
			good_marks = good_marks[hc_g$order]
			lus_smoo = lus_smoo[c(good_marks, lateral_marks),]
		}

		image(t(lus_smoo), col=mcp_heatmap_fp_shades, xaxt='n', yaxt='n')

		vals_range = range(lus_smoo)
		
		cell_x = rep(NA, time=length(cell_ord))
		names(cell_x) = names(mc@mc)
		cell_x[cell_ord] = 1:length(cell_ord)
		cl_x = tapply(cell_x, mc@mc, mean)/length(cell_ord)
		cl_x = cl_x[as.character(mc_ord)]
		mtext(mc_ord, side = 3, at=cl_x, las=2, line = 2, cex=mcp_heatmap_text_cex)
		cl_x_b = tapply(cell_x, mc@mc, max)/length(cell_ord)
		abline(v=cl_x_b, lwd=mcp_heatmap_lwd)
	} else {
		mat = log2(gene_folds[c(good_marks, lateral_marks), mc_ord])
		mat = pmax(pmin(mat,3),-3)
		if(reorder_marks) {
			hc_g = hclust(tgs_dist(mat[good_marks,]),"ward.D2")
			good_marks = good_marks[hc_g$order]
			mat = mat[c(good_marks, lateral_marks),]
		}
		image(pmax(pmin(t(mat),fold_burn), -fold_burn), col=mcp_heatmap_fp_shades, xaxt='n', yaxt='n', zlim=c(-fold_burn,fold_burn))
		mtext((1:length(mc_ord))[mc_ord], side = 3, at=seq(0,1,l=n_mc), las=2, line = 2, cex=mcp_heatmap_text_cex)
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

	if(plot_meta) {
		lower_marg=c(0,13,0,20)	
	} else {
		lower_marg=c(5,13,0,20)	
	}
	par(mar=lower_marg)
	if(plot_cells) {
		#image(as.matrix(1:length(cell_ord),nrow=1), col=mc@colors[sort(mc@mc)], xaxt='n', yaxt='n')
		if(!is.null(focus_mcs)) {
			mc_fact = as.numeric(as.factor(mc@mc[cell_ord]))
			cols = mc@colors[sort(names(mc@colors))]
			image(as.matrix(mc_fact,nrow=1), col=cols, xaxt='n', yaxt='n')
		} else {
			image(as.matrix(mc@mc[cell_ord],nrow=1), col=mc@colors, xaxt='n', yaxt='n')
		}
	} else {
		image(as.matrix(mc_ord,nrow=1), col=mc@colors, yaxt='n', xaxt='n')
	}
	if(plot_meta) {
		if(!is.null(ext_metadata)) {
			md_fact = as.character(ext_metadata[colnames(mat)])
		} else {
			md_fact = as.character(scmat@cell_metadata[colnames(mat),add_metadata])
		}
		if(sum(is.na(md_fact))>0) {
			md_fact[is.na(md_fact)] = "NA"
		}
		if(is.null(md_level_colors)) {
			md_level_colors = c("white",rep("black",times=n_md_levels))
		} else {
			md_level_colors = c("white", md_level_colors)
		}
		md_fact = as.factor(md_fact)
		n_md_levels = length(levels(md_fact))
		metadata = diag(1:n_md_levels)[md_fact,]
		message("plotting metadata on ", ncol(metadata), " rows")
		par(mar=c(5,13,0,20))
		md_hc = hclust(tgs_dist(t(metadata)),"ward.D2")
		image(metadata[,md_hc$order],xaxt='n', yaxt='n', col=md_level_colors)
		md_nms = substr(as.character(levels(md_fact)),1,20)[md_hc$order]
		mtext(md_nms,
			at=seq(0,1,length.out=ncol(metadata)),
			side=2, las=2, cex=mcp_heatmap_text_cex)
		mtext(md_nms,
			at=seq(0,1,length.out=ncol(metadata)),
			side=4, las=2, cex=mcp_heatmap_text_cex)
	}
	if(plot_cells) {
		mtext(mc_ord, side = 1, at=cl_x, las=2, line = 2, cex=mcp_heatmap_text_cex)
#		mtext(1:n_mc, side = 1, at=cl_x, las=2, line = 2, cex=mcp_heatmap_text_cex)
	} else {
		mtext((1:length(mc_ord))[mc_ord], side = 1, at=seq(0,1,l=n_mc), las=2, line = 2, cex=mcp_heatmap_text_cex)
	}
	dev.off()
	
	plot_color_bar(seq(vals_range[1], vals_range[2], length=1000), mcp_heatmap_fp_shades, gsub(".png", "_leg.png", fig_fn), show_vals_ind=c(1,500, 1000))
	
	write.table(rev(all_marks), gsub(".png", "_gene_names.txt", fig_fn), quote=F, sep="\t")
	
}

#' plot a marker heat map for a subset of metacells - selecting relevant genes for separation
#'
#' @param mc_id id of metacell object ina scdb
#' @param mat_id matrix object to us (default is the mc_id - assuming they use the same ID), if some cells in mc@mc are missing from the matrix, the function will generate an error
#' @param foc_mcs list of metacell ids to focus on
#' @param fig_fn file name for the figure (if null it will be call heat_marks in the fig directory)
#' @param n_max_mark number of markers for separation on each metacell
#' @param ddd_genes specific list of genes to add to the marker list
#' @param zero_median should median expression be used for trimming low umi counts
#' @param h height
#' @param w height

mcell_mc_plot_submc_marks = function(mc_id, mat_id, foc_mcs, fig_fn, n_max_marks=10, add_genes=NULL, zero_median=F, h = 800, w= 600)
{

	mc = scdb_mc(mc_id)
	mc = mc_restrict(mc, foc_mcs)
	lfp = log2(mc@mc_fp)
	alfp = abs(log2(mc@mc_fp))
	min_gene_fold = 1
   marks = unique(as.vector(unlist(
         apply(alfp,
            2,
            function(x)  {
               names(head(sort(-x[x>min_gene_fold]),n=n_max_marks)) })
           )))
	message("got ", length(marks))
	if(!is.null(add_genes)) {
		marks = c(add_genes, marks)
	}

	hc = hclust(tgs_dist(tgs_cor(t(lfp[marks,]))), "ward.D2")
	marks = marks[hc$order] 

	tgconfig::set_param("mcp_heatmap_height", h, "metacell")
	tgconfig::set_param("mcp_heatmap_width", w, "metacell")

	mcell_mc_plot_marks(mc_id, gset_id = mat_id, mat_id = mat_id,
								focus_mcs = foc_mcs, gene_list = marks,
								fig_fn = fig_fn, zero_median=zero_median)
}

#' Utility plot function to compare two genes based on mc_fp 
#'
#' @param mc_id metacell id in scdb
#' @param g1 gene name
#' @param g2 gene name
#' @param cex point size
#' @param text_cex text size (default 0.5)
#' @param mc_colors coloring vector per MC, default is NULL in which case we use the MC colors
#' @param fig_fn file name for the figure (optional - default is current device)
#' @param cell_md  a factor/logical vector defining which cells are in the focus. If this is unspecified all cells/metacells will be included
#' @param md_mode one of "filter", "highlight"
#' @param use_egc plot log of absolute expression instead of the log of relative expression (normalized to the median over all MCs)
#' @param md metadata to use for highlighting/filtering
#' @param md_regexp regular expression to filter/highlight metadata
#' @param md_field metadata field to use
#'
#' @export


mcell_mc_plot_gg = function(mc_id, g1, g2, 
									cex=2, text_cex=0.5, 
									mc_colors=NULL, fig_fn = NULL, 
									cell_md = NULL, md_mode="filter", 
									use_egc = F,
									main=NULL,
									cex.main=1,
									xlim=NULL, ylim=NULL,
									add_grid = F,
									e_gc_eps = 1e-5,
									md = NULL, md_regexp = NULL, md_field = NULL)
{
	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
		stop("undefined meta cell object " , mc_id)
	}
	if(is.null(mc_colors)) {
		mc_colors = mc@colors
	}
	if(use_egc) {
		lfp = log2(mc@e_gc+e_gc_eps)
	} else {
		lfp = log2(mc@mc_fp)
	}
	mc_outlines = "black"
	lwd=0.5
	mc_filt = rep(T, ncol(lfp))
	if(!is.null(md) & !is.null(md_regexp) & !is.null(md_field)) {
		cell_md = grepl(md_regexp, md[, md_field])
		names(cell_md) = rownames(md)
	}
	if(!is.null(cell_md)) {
		mc_md = 1+floor(99*tapply(cell_md[names(mc@mc)], mc@mc, mean))

		if(md_mode == "filter") {
			mc_filt = mc_md >50
		} else {
			out_shades = colorRampPalette(c("black", "darkgoldenrod", "yellow"))(100)
			mc_outlines = out_shades[mc_md]
			lwd=ifelse(mc_md > 0.25, 2, 0.5)
		}
	} 
	if(!is.null(fig_fn)) {
		png(fig_fn, w=600, h=600)
	}
	plot(lfp[g1,mc_filt], lfp[g2,mc_filt], cex=cex, xlab=g1, ylab=g2, pch=21, bg=mc_colors[mc_filt], col=mc_outlines, lwd=lwd, xlim=xlim, ylim=ylim, main=main, cex.main=cex.main)
	if(!is.na(text_cex) & text_cex != 0) {
		text(lfp[g1,mc_filt], lfp[g2,mc_filt], (1:ncol(lfp))[mc_filt], cex=text_cex)
	}
	if(add_grid) {
		grid()
	}

	if(!is.null(fig_fn)) {
		dev.off()
	}
}

#' Utility plot function to compare bulk expression of two batches/metadata factors 
#'
#' @param mat_id matrix id to use
#' @param meta_field metadata field to select cells, if not specified than a metacell must be specified
#' @param v1 first set metadata value
#' @param v2 secod set metadata value
#' @param mc_id metacell model (must be specified if meta_field is not used)
#' @param mc1 metacell id to define cells in first set
#' @param mc2 metacell id to define cells in second set
#' @param mc_to_focus metacell IDs to restrict analysis of metadata comparison
#' @param genes genes to highlight (optional)
#'
#' @export

mcell_compare_bulk = function(mat_id, 
										meta_field=NULL, v1=NULL, v2=NULL, 
										mc_id=NULL, mc1=-1, mc2=-1, 
										mc_to_focus=NULL,
										genes = NULL,  cex=0.5, text_cex=1,
										n_reg=5, fig_nm=NULL,
										top_genes=10,
										fig_h=600, fig_w=800)
{
	mat = scdb_mat(mat_id)
	if(is.null(mat)) {
		stop("undefined matrix object " , mat_id)
	}
	if(is.null(meta_field)) {
		if(is.null(mc_id)) {
			stop("must specify either meta_field or a metacell object")
		}
		if(is.null(mc1) | is.null(mc2)) {
			stop("must specify mc1 and mc2 to select cells for comparison")
		}
		mc = scdb_mc(mc_id)
		if(is.null(mc)) {
			stop("undefined meta cell object " , mc_id)
		}
		nms1 = names(mc@mc)[which(mc@mc %in% mc1)]
		nms2 = names(mc@mc)[which(mc@mc %in% mc2)]
		nms1 = intersect(nms1, colnames(mat@mat))
		nms2 = intersect(nms2, colnames(mat@mat))
		lab1 = paste("mc ", v1,sep="=")
		lab2 = paste("mc ", v2,sep="=")
		main="MC comparison"
	} else {
		md = mat@cell_metadata
		if(!meta_field %in% colnames(md)) {
			stop("Meta field ", meta_field, " does not exist")
		}
		nms1 = rownames(md)[as.character(md[,meta_field])==v1]
		nms2 = rownames(md)[as.character(md[,meta_field])==v2]
		nms1 = intersect(nms1, colnames(mat@mat))
		nms2 = intersect(nms2, colnames(mat@mat))
		lab1 = paste(meta_field, v1,sep="=")
		lab2 = paste(meta_field, v2,sep="=")
		main="batch compare on all"
		if(!is.null(mc_id)) {
			mc = scdb_mc(mc_id)
			if(is.null(mc)) {
				stop("undefined meta cell object " , mc_id)
			}
			if(is.null(mc_to_focus)) {
				stop("mc_to_focus must be specified if providing metacell object id")
			}
			cnms = names(mc@mc)[mc@mc %in% mc_to_focus]
			nms1 = intersect(nms1, cnms)
			nms2 = intersect(nms2, cnms)
			main="Batch cmp on MC subset"
		}
	}
	if(length(nms1) == 0 | length(nms2)==0) {
		stop("Got 0 cells in one of the bulks ", v1, " and ", v2)
	}
	tot1 = rowSums(mat@mat[,nms1])
	tot2 = rowSums(mat@mat[,nms2])
	N = min(sum(tot1), sum(tot2))
	tot1 = tot1*N/sum(tot1)
	tot2 = tot2*N/sum(tot2)
	
	mcell_plot_freq_compare(tot1, tot2,
						cex=cex, text_cex=text_cex,
						n_reg = n_reg, fig_nm,
						top_genes=top_genes, 
						fig_h=fig_h, fig_w=fig_w,
						lab1 = lab1, lab2 = lab2, main = main)
}

mcell_plot_freq_compare = function(tot1, tot2,
										cex=0.5, text_cex=1,
										n_reg=5, fig_nm=NULL,
										top_genes=10,
										highlight_genes = NULL,
										fig_h=600, fig_w=800,
										lab1="grp1", lab2="grp2", main="compare bulk",
										pt_col1 = "darkred", pt_col2 = "darkblue",
										show_gene_ticks = T)
{

	if(!is.null(fig_nm)) {
		png(fig_nm,w=fig_w,h=fig_h)
	}
	lr = log2(n_reg+tot2)-log2(n_reg+tot1)
	if(is.null(highlight_genes)) {
		par(mar=c(2,10,2,10))
	} else {
		par(mar=c(2,10,12,10))
	}
	xs = lr
	ys = log2(n_reg+(tot1+tot2)/2)
	max_lr = max(max(abs(lr)),1)
	plot(xs, ys, cex=cex, pch=19, col="black", ylab=NA, xlim=c(-max_lr, max_lr),
								xlab=paste(lab1, lab2,sep="/"))	
	top5 = names(tail(sort(lr),n=top_genes))
	bottom5 = names(head(sort(lr),top_genes))
	top5 = top5[order(ys[top5])]
	bottom5 = bottom5[order(ys[bottom5])]
	points(xs[bottom5], ys[bottom5], cex=cex+0.5, pch=19, col=pt_col1)
	points(xs[top5], ys[top5], cex=cex+0.5, pch=19, col=pt_col2)

	xmin = min(xs)
	xmax = max(xs)
	lab_xs = seq(xmin, xmax, length.out=length(highlight_genes))
	ymin = min(ys)
	ymax = max(ys)
	lab_ys = seq(ymin, ymax, length.out=top_genes)

	if(!is.null(highlight_genes)) {
		points(xs[highlight_genes], ys[highlight_genes], 
								cex=cex+0.5, pch=21, lwd=2,col="orange")
		mtext(highlight_genes, at=lab_xs, side=3, line=2, las=2, cex=text_cex)
		segments(x0 = lab_xs, y0=ymax+0.1*(ymax-ymin), 
							x1=xs[highlight_genes], y1=ys[highlight_genes], 
							xpd=T, lty=2)
		mtext(main, side=1)
	} else {
		mtext(main, side=3)
	}

	if(show_gene_ticks) {
		mtext(top5, at=lab_ys, side=4, line=2, las=2, cex=text_cex)
		mtext(bottom5, at=lab_ys, side=2, line=2, las=2, cex=text_cex)
		segments(x0 = -max_lr*1.1, y0=lab_ys, x1=xs[bottom5], y1=ys[bottom5], xpd=T)
		segments(x0 = max_lr*1.1, y0=lab_ys, x1=xs[top5], y1=ys[top5], xpd=T)
	}
	grid()
	if(!is.null(fig_nm)) {
		dev.off()
	}
}
