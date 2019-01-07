
#'Compute predictive value of MC cover as correlation of weighted MC averages of sc umis
#'
#'@param mc_id id of metacells to analyze
#'@param mat_id umi matrix
#'@param do_log should log transofmration be performed
#'@param focus_cells should correltion be computed to a specific set of cells
#'
mcell_wgtmc_pred_corr = function(mc_id, mat_id, graph_id, do_log = F, focus_cells=NULL)
{
	old_seed = .set_seed(get_param("mc_rseed"))

	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
		stop("MC-ERR: mc id ", mc_id, " is missing when compute mc prediction QC")
	}
	mat = scdb_mat(mat_id)
	if(is.null(mat)) {
		stop("MC-ERR: mat id ", mat_id, " is missing when computing mc prediction QC")
	}
	graph = scdb_cgraph(graph_id)
	if(is.null(graph)) {
		stop("MC-ERR: cell graph id ", graph_id, " is missing computing graph prediction QC")
	}

	umis = mat@mat

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

	mc_im = diag(max(mc@mc))[mc@mc,]
	rownames(mc_im) = names(mc@mc)

	mc_key = mc@mc[levels(graph@edges$mc2)]
	targ_mc = mc_key[graph@edges$mc2]
	neigh_im = as.matrix(table(graph@edges$mc1, targ_mc))
	neigh_im = neigh_im/rowSums(neigh_im)
	neigh_im1 = matrix(neigh_im, nrow=nrow(neigh_im))
	rownames(neigh_im1) = rownames(neigh_im)

	nfolds = 20

	all_cells = intersect(colnames(stat), names(mc@mc))
	folds = sample((1:nfolds),size=length(all_cells), rep=T)

	preds = c()

	for(fi in 1:nfolds) {
		#compute mc stats
		c_test = all_cells[folds == fi]
		c_train = all_cells[folds != fi]
		mcf_im =  mc_im[c_train,]
		mcf_im = t(t(mcf_im)/colSums(mcf_im))
		u_gm = stat[,c_train] %*% mcf_im
		pred_ig = neigh_im1[c_test,] %*% t(u_gm)
		if(is.null(preds)) {
			preds = t(pred_ig)
		} else {
			preds = cbind(preds, t(pred_ig))
		}
	}
#	save(preds, stat, file="mc_preds.Rda")
	preds = preds[,colSums(is.na(preds))==0]  #ignoring rare case of 0 connectivity within the train set
	pcor = allrow_cor(stat[,colnames(preds)], preds)
	
	.restore_seed(old_seed)

	return(pcor)	
}

