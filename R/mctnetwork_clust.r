
#' Clustring MCs based on the flow graph
#'
#' @param mct_id  mct network object
#'
#' @return A matrix on metacells and time with mean umi fraction over the gene module
#' @export
#'
mctnetwork_clust_flows = function(mct_id, K) 
{
	mct = scdb_mctnetwork(mct_id)
	if(is.null(mct)) {
		stop("cannot find mctnet object ", mct_id, " when trying to plot net flows anchors")
	}
	mc = scdb_mc(mct@mc_id)

	fls = mctnetwork_get_flow_mat(mct,-1)
	fls3 = fls %*% fls %*% fls
	cr3 = tgs_cor(rbind(fls3, t(fls3)))
	diag(cr3) = 0
	cr3_2 = tgs_cor(cr3)
	diag(cr3_2) = 0

	hc = hclust(tgs_dist(t(cr3_2)), "ward.D")

	clf = cutree(hc, K)
	return(list(cmat = cr3_2, clust = clf, hc=hc))
}

mctnetwork_plot_cormat = function(mct_id, mat_id, gset_id, K, fn, w=1600, h=1600, text_cex=1, cls_pre_order=NULL, marks=NULL)
{
	mct = scdb_mctnetwork(mct_id)
	if(is.null(mct)) {
		stop("cannot find mctnet object ", mct_id, " when trying to plot net flows anchors")
	}
	clst_flows = mctnetwork_clust_flows(mct_id, K)

	cmat = clst_flows$cmat
	fclst = clst_flows$clust

	mc_ord = clst_flows$hc$order
	if(!is.null(cls_pre_order)) {
		cls_rank = 1:K
		names(cls_rank) = as.character(1:K)
		cls_rank[as.character(cls_pre_order)] = 1:K
		mc_mean_age = apply(mct@mc_t, 1, function(x) {return(mean(x*(1:length(x))/sum(x))) })
		mc_ord = order(-cls_rank[fclst]*1000+mc_mean_age)
	}

	mc = scdb_mc(mct@mc_id)

	shades = colorRampPalette(c("darkblue", "blue","white", "red", "yellow"))(1000)
	png(fn, w=w, h=h)
	n_mc = nrow(cmat)

	layout(matrix(c(1,2),nrow=2),heights=c(nrow(cmat)*9+50, 300))
	par(mar=c(0,5,4,5))
	image(cmat[mc_ord, mc_ord], zlim=c(-1,1), col=shades, yaxt='n', xaxt='n')
	N = length(mc_ord)

	mc_x = 1:length(mc_ord)
	names(mc_x) = 1:N
	mc_x[mc_ord] = 1:N
	for(i in 1:K) {
		abline(h=max(-0.5+mc_x[fclst == i])/(N-1))
		abline(v=max(-0.5+mc_x[fclst == i])/(N-1))
	}
	cl_x = tapply(mc_x, fclst, mean)
	cl_max = tapply(mc_x, fclst, max)
	
	mtext(1:K, side = 3, at=cl_x/N, las=1, cex=1.5)

	mtext((1:length(mc_ord))[mc_ord], side = 2, at=seq(0,1,l=n_mc), las=2, line = 2, cex=text_cex)
	mtext((1:length(mc_ord))[mc_ord], side = 4, at=seq(0,1,l=n_mc), las=2, line = 2, cex=text_cex)
	par(mar=c(3,5,0,5))
	image(as.matrix(mc_ord,nrow=1), col=mc@colors, yaxt='n', xaxt='n')
	mtext((1:length(mc_ord))[mc_ord], side = 1, at=seq(0,1,l=n_mc), las=2, line = 2, cex=text_cex)
	dev.off()

	if(!is.null(marks)) {
		egc = log2(mc@e_gc+1e-5)
		dir = sub(".png","_marks",fn)
		if(!dir.exists(dir)) {
			dir.create(dir)
		}
		for(g in marks) {
			png(sprintf("%s/%s.png", dir, g),w=1000,h=300)
			plot(1:length(mc_ord), egc[g, mc_ord], pch=19, col=mc@colors[mc_ord], ylab=g,xaxt='n')
			mtext(1:K, at=cl_x,side=1, las=2)
			
			abline(v=cl_max+0.5)
			grid()
			dev.off()
		}
	}
	mark_fn = sub("png", "mat.png", fn)
	mcell_mc_plot_marks(mct@mc_id, gset_id=gset_id, mat_id = mat_id, fig_fn=mark_fn, mc_ord = mc_ord, plot_cells=T)
	return(mc_ord)
}	
