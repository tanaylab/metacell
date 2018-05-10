#' Simple screen for outlier cells in a metacell cover, finding genes with overly high expression given their metacell mean
#'
#'
#' @param new_mc_id id of new metacell to create
#' @param mc_id ID of metacell in scdb
#' @param mat_id matrix object to use (should be compatible with mc_id)
#' @param T_lfc maximal log2 fold change over the expected to be considered as an outlier.
#'
#' @export
#'
#' @import dbscan
#' @import entropy
#'
mcell_mc_screen_outliers_1gene_fold = function(new_mc_id, mc_id, mat_id, T_lfc)
{
	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
		stop("cannot find mc ", mc_id, " when attempting simple outlier detection")
	}
	mat = scdb_mat(mat_id)
	if(is.null(mat)) {
		stop("cannot find mat ", mat_id, " when attempting simple outlier detection")
	}

	lfc_gi = mc_compute_outlier_fc(mc, mat)
	maxlfc_g = apply(lfc_gi, 1, max)
	maxlfc_i = apply(lfc_gi, 2, max)

	if(sum(maxlfc_i > T_lfc) > 0) {
		new_outliers = names(which(maxlfc_i > T_lfc))
		mc@outliers = c(mc@outliers, new_outliers)
		mc@mc = mc@mc[setdiff(names(mc@mc),new_outliers)]
		scdb_add_mc(new_mc_id, tgMCCov(mc@mc, mc@outliers, mat))
	}
}

#' Compute log fold change expression of each gene given its regularized metacell expression
#'
#' Each metacell is defiend by a total set of umis, defining an expected number of umi pere cell per gene. After regularization (1 umi per cell), the fold enrichment of the gene in each cell is computed naively. Note that the MC inferred value is currently based on linear mean, whic his highly inadequate when an outlier is present.
#'
#' @param mc metacell object (not id)
#' @param mat matrix object (not id)
#'
#' @export
#'
mc_compute_outlier_fc = function(mc, mat)
{
	min_outlier_u = 6

	ismc_ci = diag(max(mc@mc))[,mc@mc]
	
	#we compute the expected number of umi per gene per cell given clust
	u_gi = mat@mat[,names(mc@mc)]

	u_i = colSums(u_gi)

	#saveing computation on boring genes
	ishigh_g = rep(F, times=nrow(u_gi))
	for(i in 0:floor(nrow(u_gi)/1000)) {
		fr = i * 1000 + 1
		to = min(fr+1000, nrow(u_gi))
		ishigh_g[fr:to] = apply(u_gi[fr:to,],1,max) >= min_outlier_u
	}
	u_gi = u_gi[ishigh_g,]

	u_gc = u_gi %*% t(ismc_ci)
	u_c = colSums(u_gc)

	p_gc = t(t(u_gc) / u_c)
	p_gi = p_gc %*% ismc_ci
	exp_gi = t(t(p_gi) * u_i)

	k_ureg = 1
	lfc_gi = log2((k_ureg+u_gi)/(k_ureg+exp_gi))

	return(lfc_gi)
}

#' Plot and outlier heat map. 
#'
#' This select genes that are expressed in an interesting way in one of the outliers, and visualize their combinatorics as a heatmap.
#'
#' @param mc metacell object (not id)
#' @param mat matrix object (not id)
#'
#' @export
mcell_plot_outlier_heatmap = function(mc_id, mat_id, T_lfc)
{
	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
		stop("cannot find mc ", mc_id, " when attempting simple outlier detection")
	}
	mat = scdb_mat(mat_id)
	if(is.null(mat)) {
		stop("cannot find mat ", mat_id, " when attempting simple outlier detection")
	}
	lfc_gi = mc_compute_outlier_fc(mc, mat)
	u_gi = mat@mat[,names(mc@mc)]
	maxlfc_g = apply(lfc_gi, 1, max)
	maxlfc_i = apply(lfc_gi, 2, max)

	out_g_nms = names(which(maxlfc_g > T_lfc))
	out_i_nms = names(which(maxlfc_i > T_lfc))

	if(sum(maxlfc_g > T_lfc) > 1 & sum(maxlfc_i > T_lfc) > 2) {
		outu_gi = log2(1 + as.matrix(u_gi[out_g_nms, out_i_nms]))

	#reporting the outlier gene / cell matrix

		hc1 = hclust(dist(cor(outu_gi)), "ward.D2")
		hc2 = hclust(dist(cor(t(outu_gi))), "ward.D2")
		fig_nm = scfigs_fn(mc_id, "outlier")
	
		h_mat = min(300+length(out_g_nms)*20,2000)
		png(fig_nm, w=min(300+20*length(out_i_nms),3000), h=h_mat+100)

		layout(matrix(c(1,2),nrow=2),heights=c(h_mat, 100))
		top_marg=c(0,13,5,20)
		par(mar=top_marg)
		shades = colorRampPalette(c("white", "blue", "red", "yellow", "black"))(1000)
		c_ord = order(mc@mc[out_i_nms])
		image(t(outu_gi[hc2$order, c_ord]), col=shades, xaxt='n', yaxt='n')
		mtext(rownames(outu_gi)[hc2$order], at=seq(0,1,l=length(hc2$order)), las=2, side=2, cex=0.8)
		mtext(rownames(outu_gi)[hc2$order], at=seq(0,1,l=length(hc2$order)), las=2, side=4, cex=0.8)

		lower_marg=c(5,13,0,20)
		par(mar=lower_marg)
		mc_cols = mc@colors[mc@mc[out_i_nms[c_ord]]]
		image(as.matrix(1:length(out_i_nms),nrow=1), col=mc_cols, xaxt='n', yaxt='n')
		dev.off()
	}
}


#' Split and filter metacells using dbscan and outlier gene detection
#'
#' Running intra-metacell clustering on each metacell, identifying splits if existing. Also remove cells with outlier expression from each metacell. Note that this can in many case eliminate doublets as outliers - but may be insffucieint as doublets may constitute homogeneous metacells in large datasets.
#'
#' @param new_mc_id id of metacells to create
#' @param mc_id  id of source metacell object
#' @param mat_id  id of umi matrix
#' @param plot_mats if this is true, a heatmap will be generated for each metacell
#' @param dirichlet if this is true, the distance metric to be used for splitting is dirichlet on the downsampled umi counts. Otherwise (and by default), the correlation between log(1+umi) of the downsampled matrix will be used. Note that downsampling is done given the miminum depth of each metacell, seperately.
#'
#' @export
mcell_mc_split_filt = function(new_mc_id, mc_id, mat_id, T_lfc, plot_mats=T, dirichlet=F)
{
	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
		stop("cannot find mc ", mc_id, " when attempting simple outlier detection")
	}
	mat = scdb_mat(mat_id)
	if(is.null(mat)) {
		stop("cannot find mat ", mat_id, " when attempting simple outlier detection")
	}

	if(plot_mats) {
	  plot_dir = scfigs_dir(mc_id, "outliers")
	  file.remove(Sys.glob(sprintf("%s/*", plot_dir)))
	}	

	split_dbscan = function(nms) {
		id = mc@mc[nms[1]]
		mc_mat = mat@mat[,nms]
		mc_mat = scm_downsamp(mc_mat, min(colSums(mc_mat)))
		e = rowMeans(mc_mat)
		v = apply(mc_mat, 1, var)
		gs = (v+0.1)/(e+0.1) > 1.2 & apply(mc_mat, 1, max)>2
		if(sum(gs) < 2 | sum(colSums(mc_mat[gs,] > 3)) < 20) {
			message("no genes for ", id)
			clst = rep(1, times=length(nms))
			names(clst) = nms
			return(clst)
		}
		mc_mat = mc_mat[gs,]
		null_nms = names(which(colSums(mc_mat) <= 3 | apply(mc_mat,2,var)==0))
		filt_nms = setdiff(nms, null_nms)
		if(length(null_nms) > 0) {
			message("got ", length(null_nms), " cells w.o enough umis")
		}
		mc_mat = mc_mat[, filt_nms]
		
		if(dirichlet) {
			reg = 0.5
			L = length(nms)
			dst = matrix(NA, nrow=L, ncol=L)
			for(i in 1:(L-1)) {
				for(j in (i+1):L) {
					dst[i,j] = KL.Dirichlet(mc_mat[,i], mc_mat[,j], reg, reg)
				}
			}
			dst[is.na(dst)] = t(dst)[is.na(dst)]
			diag(dst) = 0
		} else {
			dst = 1-tgs_cor(log(1+7*as.matrix(mc_mat)))
		}
		clst = dbscan(as.dist(dst), eps=quantile(dst,0.1), minPts=5)
		clst = clst$cluster
		names(clst) = filt_nms
		clst[null_nms] = 0
		if(is.na(sum(clst<2))) {
			message("NA after dbscan, clst id ", id)
			clst = rep(1, times=length(nms))
			names(clst) = nms
			return(clst)
		}
#		if(sum(clst!=1) > 0) {
#			message("breaking ", id, " : ", 
#					paste(as.numeric(table(clst)), collapse=" "))
#		} else {
#			message("homogeneous ", id)
#		}
		clst[is.na(clst)] = 0
		return(clst)
	}
	mc_cores = get_param("mc_cores")
	nms_mc = split(names(mc@mc), mc@mc)
	message("starting split outlier dbscan")
#	all_clst = lapply(nms_mc, split_dbscan)
	all_clst = mclapply(nms_mc, split_dbscan, mc.cores=mc_cores)

	new_mc = rep(NA, length(mc@mc))
	names(new_mc) = names(mc@mc)
	next_mcid = 1
	good_cells = names(new_mc)
	clust_outliers = c()
	for(cid in 1:length(all_clst)) {
#		message("opening cid ", cid)
		clst = all_clst[[cid]]
#		message("opening clst length ", length(clst))
		#adding cluster 1 and all dbscan outliers that were not filtered parametrically to a new MC
		new_mc[intersect(good_cells, names(which(clst < 2)))] = next_mcid
		next_mcid = next_mcid + 1
		if(length(clst) > 1 & max(clst) > 1) {
		  for(sub_i in 2:max(clst)) {
				message("splitting metacell ", cid)
				new_mc[intersect(good_cells, names(which(clst == sub_i)))] = next_mcid
				clust_outliers = c(clust_outliers, 
											names(which(clst == sub_i)))
				next_mcid = next_mcid + 1
		  }  
		}
	}
	lfc_gi = mc_compute_outlier_fc(mc, mat)
	if(plot_mats) {
		for(clst in all_clst) {
			mc_plot_submc(mc, clst, mat, plot_dir, lfc_gi[,names(clst)])
		}
	}
	maxlfc_i = apply(lfc_gi, 2, max)

	new_outliers = setdiff(names(which(maxlfc_i > T_lfc)), clust_outliers)
	new_outliers = c(new_outliers, names(which(is.na(new_mc))))
	new_outliers = c(mc@outliers, new_outliers)
	new_mc = new_mc[setdiff(names(new_mc),new_outliers)]
	scdb_add_mc(new_mc_id, tgMCCov(new_mc, new_outliers, mat))
}

mc_plot_submc = function(mc, clst, mat, base_dir, lfc_gi)
{
	nms = names(clst)
	mc_id = mc@mc[nms[1]]
	lfp = log2(mc@mc_fp[,mc_id])

	mc_genes = names(tail(sort(lfp),5))

	us = mat@mat[,nms]
	top_genes = union(names(which(apply(us,1,max)>3)),mc_genes)
	us = us[top_genes,]
	sub_fp = t(apply(us,1, function(x) tapply(x,clst,mean)))

	sub_fp = log2((0.1+sub_fp)/(0.1+rowMeans(sub_fp)))
	max_fc0 = sub_fp[,1]
	max_fc0 = max_fc0[max_fc0 > 0.3]
	max_fc0 = sort(max_fc0)
	noise_genes = names(tail(max_fc0, 10))
	fc = apply(abs(sub_fp), 1, max)
	max_fc = fc[fc > 0.5]
	max_fc = sort(max_fc)
	out_genes = c(noise_genes, names(tail(max_fc, 5)))

	genes = unique(c(out_genes, mc_genes))

	lfc_genes = names(which(apply(lfc_gi,1,max) > 3))
	lfc_max = apply(lfc_gi,2, max)
	genes = c(lfc_genes,setdiff(genes,lfc_genes))

	us = us[genes,order(clst + lfc_max/100)]

	png(sprintf("%s/mc%d_out.png", base_dir, mc_id), 
											w=400+ncol(us)*5, h=400+nrow(us)*20)
	top_marg=c(5,13,5,20)
	par(mar=top_marg)
	shades = colorRampPalette(c("white", "blue", "red", "yellow", "black"))(1000)
	desc = sprintf("noise %d, clsts %s ", sum(clst==0), 
						paste(as.numeric(table(clst[clst>0])),collapse=" "))
	image(log2(1+t(as.matrix(us))), col=shades, xaxt='n', yaxt='n', main=desc)
	abline(v=(-0.5+sum(clst==0))/(length(clst)-1), lwd=3)
	abline(h=(-0.5+length(lfc_genes))/(nrow(us)-1), lwd=3)
	mtext(rownames(us), at=seq(0,1,l=nrow(us)), las=2, side=2, cex=0.8)
	mtext(rownames(us), at=seq(0,1,l=nrow(us)), las=2, side=4, cex=0.8)
	dev.off()
}
