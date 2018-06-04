#`compute MC cgraph edge density
#'
#' for each MC we compute mean and variance of the intra-MC edge density
#'
#'@param mc_id MC cover to analyze
#'@param cgraph_id which cgraph to study
#'
#'
mcell_mc_edge_density = function(mc_id, graph_id, mc2d_id = NA)
{
	graph = scdb_cgraph(graph_id)
	if(is.null(graph)) {
		stop("MC-ERR: cell graph id ", graph_id, " is missing when mc 2d projection")
	}
	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
		stop("MC-ERR: mc id ", mc_id, " is missing when running add_mc_from_graph")
	}

	c_homog = mcell_mc_cell_homogeneity(mc_id, graph_id, ignore_mis=T)

	cnms = rownames(c_homog)
	c_homog = c_homog[,2]/rowSums(c_homog)
	homog_e = tapply(c_homog, mc@mc[cnms], mean)
	homog_sd = tapply(c_homog, mc@mc[cnms], sd)

	fig_nm = scfigs_fn(mc_id, sprintf("cgraph_%s_homog_e_sd", graph_id))
	png(fig_nm,w=800,h=800)
	plot(homog_e, homog_sd, pch=19, col=mc@colors, cex=3, xlab="mean homogeneity", ylab="sd homogeneity")
	text(homog_e, homog_sd, 1:max(mc@mc), cex=0.8)
	dev.off()
	fig_nm = scfigs_fn(mc_id, sprintf("cgraph_%s_mc_sz_homog", graph_id))
	png(fig_nm,w=800,h=800)
	plot(tabulate(mc@mc), homog_e, pch=19, col=mc@colors, cex=3, xlab="cell number", ylab="mean homogeneity")
	text(tabulate(mc@mc), homog_e, 1:max(mc@mc), cex=0.8)
	dev.off()
	
}

#'Compute predictive value of MC cover as correlation of MC averages of sc umis
#'
#'@param mat_id id of matrix to analyze
#'@param graph_id id of cgraph to use for computing neighborhood per cell
#'@param do_log should log transofmration be performed
#'@param focus_cells should correltion be computed to a specific set of cells
#'
mcell_graph_pred_corr = function(mat_id, graph_id, do_log = F, focus_cells=NULL, filt_genes=NULL)
{
	graph = scdb_cgraph(graph_id)
	if(is.null(graph)) {
		stop("MC-ERR: cell graph id ", graph_id, " is missing computing graph prediction QC")
	}
	mat = scdb_mat(mat_id)
	if(is.null(mat)) {
		stop("MC-ERR: mat id ", mat_id, " is missing when computing graph prediction QC")
	}

	umis = mat@mat
	if(!is.null(filt_genes)) {
		umis = umis[setdiff(rownames(umis),filt_genes),]
	}

	csize = colSums(umis)
	n_downsamp = quantile(csize,0.05)
	umis_ds = scm_downsamp(umis, n_downsamp)

	f = rowSums(umis) > 100	#just to save comp time
	umis = umis[f,]
	umis_ds = umis_ds[f,]
	if(do_log) {
		stat = Matrix(log2(1+ 7*umis_ds))
	} else {
		stat = umis_ds
	}

	edges = graph@edges

	N = length(levels(edges$mc1))
	cell_neig = sparseMatrix(i=as.integer(edges$mc1), 
							j=as.integer(edges$mc2), x=1, dims=c(N,N))
	colnames(cell_neig) = as.character(levels(edges$mc1))
	diag(cell_neig) = 0
	joint_cnms = intersect(as.character(levels(edges$mc1)), colnames(stat))

	if(!is.null(focus_cells)) {
		joint_cnms = intersect(focus_cells, joint_cnms)
	}

	cell_neig = cell_neig[,joint_cnms]

	cell_neig = cell_neig/rowSums(cell_neig)

	stat = stat[,joint_cnms]

	stat_neig = t(cell_neig %*% t(stat))

	colnames(stat_neig) = as.character(levels(edges$mc1))

	corstat = stat; corstat_neig = stat_neig[,joint_cnms]
	save(corstat_neig, corstat, file="corr_preds.Rda")
	pred = allrow_cor(stat, stat_neig[,joint_cnms])
	return(pred)	
}

#'Compute predictive value of MC cover as correlation of MC averages of sc umis
#'
#'@param mc_id id of metacells to analyze
#'@param mat_id umi matrix
#'@param do_log should log transofmration be performed
#'@param focus_cells should correltion be computed to a specific set of cells
#'
mcell_mc_pred_corr = function(mc_id, mat_id, do_log = F, focus_cells=NULL, filt_genes=NULL)
{
	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
		stop("MC-ERR: mc id ", mc_id, " is missing when compute mc prediction QC")
	}
	mat = scdb_mat(mat_id)
	if(is.null(mat)) {
		stop("MC-ERR: mat id ", mat_id, " is missing when computing mc prediction QC")
	}

	umis = mat@mat
	if(!is.null(filt_genes)) {
		umis = umis[setdiff(rownames(umis),filt_genes),]
	}

	csize = colSums(umis)
	n_downsamp = quantile(csize,0.05)
	umis_ds = scm_downsamp(umis, n_downsamp)

	f = rowSums(umis_ds) > 100	#just to save comp time
	umis = umis[f,]
	umis_ds = umis_ds[f,]
	if(do_log) {
		stat = Matrix(log2(1+ 7*umis_ds))
	} else {
		stat = umis_ds
	}
	l = lapply(tabulate(mc@mc), function(n) matrix(rep(1/n, n*n), nrow=n))
 	cell_neig = .bdiag(l)
	colnames(cell_neig) = names(mc@mc)[order(mc@mc)]
	rownames(cell_neig) = names(mc@mc)[order(mc@mc)]
	diag(cell_neig) = 0

	joint_cnms = intersect(names(mc@mc), colnames(stat))
	if(!is.null(focus_cells)) {
		joint_cnms = intersect(focus_cells, joint_cnms)
	}
	cell_neig = cell_neig[,joint_cnms]

	cell_neig = cell_neig/rowSums(cell_neig)

	stat = stat[,joint_cnms]

	stat_neig = t(cell_neig %*% t(stat))

	colnames(stat_neig) = rownames(cell_neig)

	mcstat = stat; mcstat_neig = stat_neig[,joint_cnms]
	save(mcstat_neig, mcstat, file="mc_preds.Rda")
	pred = allrow_cor(stat, stat_neig[,joint_cnms])
	return(pred)	
}

