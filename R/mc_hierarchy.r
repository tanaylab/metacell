#' comput mc hierarchucal clustering using the normalized confusion matrix
#'
#' @param mc_id metacell id
#' @param graph_id cgraph to define edges and impose them on the metacells
#' @param confu a confusion matrix derived from e.g. mcell_mc_confusion_matm or mcell_mc_coclust_confusion_mat
#'
#' @export

mcell_mc_hclust_confu = function(mc_id, graph_id, confu=NULL)
{
	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
		stop("undefined meta cell object " , mc_id)
	}
	if(is.null(graph_id)) {
		if(is.null(confu)) {
			stop("MC-ERR, missing graph_id AND confu object when trying to hclust mc by confusion matrix")
		}
		if(nrow(confu) != ncol(mc@mc_fp) | ncol(confu) != nrow(confu)) {
			stop("MC-ERR, bad ocnfusion matrix dimension when trying to hcluster mc by confusion matrix")
		}
	} else {
		cgraph = scdb_cgraph(graph_id)
		if(is.null(cgraph)) {
			stop("undefined cgraph object when trying to plot confusion, id " , graph_id)
		}
	
		max_deg = nrow(cgraph@edges)
		confu = mcell_mc_confusion_mat(mc_id, graph_id, max_deg, ignore_mismatch=T)
	}
	r_confu = rowSums(confu)
	c_confu = colSums(confu)
	norm = r_confu %*% t(c_confu)
	confu_n = confu/norm

	confu_nodiag = confu_n
	diag(confu_nodiag) = 0
	confu_n = pmin(confu_n, max(confu_nodiag))
	confu_n = pmin(confu_n, quantile(confu_n, 1-3/nrow(confu_n)))
	epsilon = quantile(confu_n[confu_n!=0],0.02)
	hc = hclust(as.dist(-log10(epsilon+confu_n)),"average")
	return(hc)
}

#' identify super structure in an mc cover, based on hcluster of the confusion matrix
#'
#' @param mc_id id of metacell object ina scdb
#' @param mc_hc hclust object onthe metacells (e.g. derive from mcell_mc_hclust_confu)
#' @param T_gap the minimal branch length for defining a supper metacell structure
#'
#' @export

mcell_mc_hierarchy = function(mc_id, mc_hc, T_gap)
{
	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
		stop("undefined meta cell object " , mc_id)
	}
	
	mc_ord = rep(0, length(mc_hc$order))
	mc_ord[mc_hc$order] = 1:length(mc_hc$order)

	n_mc = ncol(mc@mc_fp)

	lfp = log2(mc@mc_fp)

	parent = rep(-1, n_mc-1)
	cells = rep(0, times=n_mc-1)
	mcs = list()
	gaps = rep(0, times=n_mc-1)
	for(i in 1:nrow(mc_hc$merge)) {
		left = mc_hc$merge[i,1]
		right = mc_hc$merge[i,2]
		if(left < 0) {
			cs = names(which(mc@mc == -left))
			cell_left = length(cs)
			mc_left = -left
		} else {
			cell_left = cells[left]
			mc_left = mcs[left]
			parent[left] = i
			gaps[left] = mc_hc$height[i] - mc_hc$height[left]
		}
		if(right < 0) {
			cs = names(which(mc@mc == -right))
			cell_right= length(cs)
			mc_right = -right
		} else {
			cell_right= cells[right]
			mc_right = mcs[right]
			parent[right] = i
			gaps[right] = mc_hc$height[i] - mc_hc$height[right]
		}
#		message("in node ", i, " l ", left, " r ", right)
		mcs[[i]] = c(unlist(mc_right), unlist(mc_left))
		cells[i] = cell_right + cell_left
	}
	n_min_outcells = 300
	hits = list()
#	covered = rep(0, times = length(gaps))
	hit_i = 1
	sup_x = c()
	for(i in which(gaps > T_gap)) {
		j = parent[i]
		mincells = cells[i] + n_min_outcells
		while(j != -1 & cells[j] < mincells) {
			j = parent[j]	
		}
		if(j != -1) { 
#			message("node ", i, " sup ", j)
#& covered[j] == 0
			N = cells[j]
			n = cells[i]
			mcs_in = mcs[[i]]
			mcs_out = setdiff(mcs[[j]], mcs[[i]])
			lfp_avg_in = apply(lfp[,mcs_in], 1, mean)
			lfp_min_in = apply(lfp[,mcs_in], 1, min)
			lfp_max_in = apply(lfp[,mcs_in], 1, max)
			if(length(mcs_out) > 1) {
				lfp_max_out= apply(lfp[,mcs_out], 1, max)
				lfp_min_out= apply(lfp[,mcs_out], 1, min)
				lfp_avg_out= apply(lfp[,mcs_out], 1, mean)
			} else {
				if(length(mcs_out) == 0) {
					message("zero length mcs out??")
					stop("boom")
				}
				lfp_max_out= lfp[,mcs_out]
				lfp_min_out= lfp[,mcs_out]
				lfp_avg_out= lfp[,mcs_out]
			}

			e_marks = tail(sort(lfp_avg_in),20)
			sep_marks = tail(sort(lfp_min_in),20)
			marks_gap = tail(sort(lfp_avg_in - lfp_avg_out),20)
			marks_gap_anti = head(sort(lfp_avg_in - lfp_avg_out),20)
			
			x_ord = mean(mc_ord[mcs[[i]]])
			sup_x = c(sup_x, x_ord)
			hits[[hit_i]] = list(marks = e_marks, 
								min_marks = sep_marks,
								marks_gap = marks_gap, 
								marks_gap_anti=marks_gap_anti, 
								mcs = mcs[[i]], 
								x_ord = x_ord, sup_mcs = mcs[[j]])
			hit_i = hit_i + 1
		}
	}

	hits = hits[order(sup_x)]
	return(hits)
}

#' plot super strucutre: super clust mc footprint, and selected genes
#'
#' @param mc_id id of metacell object
#' @param graph_id id of graph for (Re-) constructing confusion matrix
#' @param mc_order the mc ordering (e.g., hc$order using the output of mcell_mc_hclust_confu)
#' @param sup_mcs the list you get from mcell_mc_hierarchy (for now)
#' @param width width of figure in pixels
#' @param height height of figure in pixels
#' @param fig_fn figure name (NULL will create a figure named [mc_id]_supmc_confu in the figure directory)
#' @param min_nmc minimal number of mc in supmc set, smaller gruops will not be plotted
#' @param shades heatmap color palette 
#' @param plot_grid plot vertical grid in the heatmap
#' @param show_mc_ids plot mc ids below the heatmap
#'
#' @export

mcell_mc_plot_hierarchy = function(mc_id, graph_id, mc_order, sup_mcs, width, height, fig_fn=NULL, min_nmc=2, shades = colorRampPalette(c("white", "pink", "red", "black", "brown", "orange")), plot_grid=T, show_mc_ids=F)
{
	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
		stop("undefined meta cell object " , mc_id)
	}
	n_mc = ncol(mc@mc_fp)
	colors = mc@colors

	confu = mcell_mc_confusion_mat(mc_id, graph_id, 1000, ignore_mismatch=T)
	r_confu = rowSums(confu)
	c_confu = colSums(confu)
	norm = r_confu %*% t(c_confu)
	confu_n = confu/norm

	confu_n = confu_n[mc_order, mc_order]
	colors = colors[mc_order]
	colnames(confu_n) = (1:ncol(confu_n))[mc_order]
	rownames(confu_n) = (1:ncol(confu_n))[mc_order]
	colors[is.na(colors)] = "gray"

	fps = as.matrix(do.call(cbind, lapply(sup_mcs, 
					function(x) { fp = rep(0, n_mc); 
										fp[unlist(x$sup_mcs)]=1; 
										fp[unlist(x$mcs)]=2; 
										fp })))


	if(is.null(fig_fn)) {
		fig_fn = scfigs_fn(mc_id, "supmc_confu")
	}
	marks = unlist(lapply(sup_mcs, function(x) {
						m = x$marks[x$marks > 0.4]
						return(paste(names(tail(m,5)),collapse=","))
					}
				))

	gap_marks = unlist(lapply(sup_mcs, function(x) {
						m = x$marks_gap[x$marks_gap > 0.4];
						s = rev(names(tail(m,5)));
						return(substr(paste(s, collapse=","),1,30))
					}
					))
	gap_amarks = unlist(lapply(sup_mcs, function(x) {
						m = x$marks_gap_anti[x$marks_gap_anti < -0.4];
						s = rev(names(head(m,5)));
						return(substr(paste(s, collapse=","),1,30))
					}
					))
	gap_marks = paste(gap_marks,gap_amarks,sep=" | ")
	marks = paste(marks, 1:length(marks), sep=" :")
	vert = F

	f_sup = colSums(fps>1)>min_nmc

	grids = -1+apply(fps[mc_order,f_sup],2,function(x) min(which(x>1)))
	grids = c(grids, apply(fps[mc_order,f_sup],2,function(x) max(which(x>1))))

	marks = marks[f_sup]
	gap_marks = gap_marks[f_sup]

	.plot_start(fig_fn, w=width, h=height)

	layout(matrix(c(1,2,3), ncol=1), heights=c(10,10,0.5))

	par(mar=c(0,30,2,40))
	image(fps[mc_order, f_sup], col=c("white", "lightgray", "blue"), xaxt='n', yaxt='n')
	mtext(marks, side = 2, at=seq(0,1,length.out=length(marks)), las=2, cex=1)
	mtext(gap_marks, side = 4, at=seq(0,1,len=length(marks)), las=2, cex=1)
	if (plot_grid) {
		abline(v=(-0.5+grids)/(n_mc-1), lwd=0.5)
	}

	par(mar=c(0,30,0,40))
	confu_nodiag = confu_n
	diag(confu_nodiag) = 0
	confu_n = pmin(confu_n, max(confu_nodiag))
	confu_n = pmin(confu_n, quantile(confu_n, 1-3/nrow(confu_n)))
	image(confu_n,col=shades(1000),xaxt='n', yaxt='n')
	if (plot_grid) {
		abline(v=(-0.5+grids)/(n_mc-1), lwd=0.5)
	}
	par(mar=c(5,30,0,40))
	image(as.matrix(1:n_mc,nrow=1), col=colors, yaxt='n', xaxt='n')

	if (show_mc_ids) {
		mtext(colnames(confu_n), side=1, at=seq(0, 1, len=ncol(confu_n)), las=2, line=1, cex=0.7)
	}
	dev.off()
}
