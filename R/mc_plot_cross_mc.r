#' Plotting a matrix of co-occurences between two metacell covers of the same dataset
#'
#' @param mc_id1 first mc object to use
#' @param mc_id2 second mc object to use
#' @param ord_1 optional ordering of first MC set (default is original order)
#' @param ord_2 optional ordering of second MC set (Default is based on max intersection with the first set)
#' @param not_log use this to avoid using log scale on the colors - making the best match MC more pronanouced
#'
#' @export
mcell_plot_cross_mc = function(mc_id1, mc_id2, ord_1=NULL, ord_2=NULL, not_log=F)
{
	mc1 = scdb_mc(mc_id1)
	mc2 = scdb_mc(mc_id2)
	if(is.null(mc1) | is.null(mc2)) {
		stop("missing metacell ids = ", mc_id1, " or ", mc_id2)
	}
	n1 = max(mc1@mc)
	n2 = max(mc2@mc)
	joint_c = intersect(names(mc1@mc), names(mc2@mc))
	i_mc1 = mc1@mc[joint_c]
	i_mc2 = mc2@mc[joint_c]
	all=union(names(i_mc1), names(i_mc2))
	i_mc1[setdiff(all,names(i_mc1))] = 0
	i_mc2[setdiff(all,names(i_mc2))] = 0
	cross = t(matrix(tabulate(i_mc1 * (n2+1) + i_mc2 + 1, nbins=(1+n1)*(1+n2)), nrow=(n2+1)))

#	hc = hclust(dist(cor(cross)), "ward.D2")
	if(is.null(ord_1)) {
		ord_1 = 1:nrow(cross)
	} else {	#adding the 0 row
		ord_1 = c(1, ord_1+1)
	}
	if(is.null(ord_2)) {
		ord_2 = order(apply(cross[ord_1,],2, which.max))
	} else { #adding the 0 column
		ord_2 = c(1, ord_2+1)
	}

	fig_nm = scfigs_fn(mc_id1, sprintf("_cross_%s", mc_id2))
	shades = colorRampPalette(c("white", "pink", "red", "black", "brown", "orange"))(1000)
	width = min(n1*50+400,3000)
	height = min(n2*50+400,3000)
	png(fig_nm, width = width, height = height)
	layout(matrix(c(1,4,2,3),nrow=2),
					widths=c(150, width-150),
					heights=c(height-100, 100))
	top_marg=c(0,8,4,0)
	par(mar=top_marg)
	image(t(as.matrix(1:(n2+1),nrow=1)), col=(c("white",mc2@colors))[ord_2], yaxt='n', xaxt='n')
	mtext((0:n2)[ord_2], at=seq(0,1,l=n2+1), las=2, side=2, cex=1.5)
	top_marg=c(0,0,4,4)
	par(mar=top_marg)
	if(not_log) {
		image(log2(1+cross)[ord_1,ord_2], col=shades, xaxt='n', yaxt='n')
	} else {
		image(log2(1+cross)[ord_1,ord_2], col=shades, xaxt='n', yaxt='n')
	}
	lower_marg=c(4,0,0,4)
	par(mar=lower_marg)
	image(as.matrix(1:(n1+1),ncol=1), col=c("white", mc1@colors)[ord_1], yaxt='n', xaxt='n')
	mtext((0:n1), at=seq(0,1,l=n1+1), las=2, side=1, cex=1.5)
	dev.off()
}
