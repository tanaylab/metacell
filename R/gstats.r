
#' mcell_gene_stat
#'
#' @param mat_id
#'
#' @export

mcell_add_gene_stat = function(mat_id, gstat_id, force=F)
{
	if(!scdb_is_valid()) {
		stop("MCERR - scdb is not initialized, cannot import")
	}
	if(!force & !is.null(scdb_gstat(gstat_id))) {
		return(TRUE)
	}
	if(is.null(scdb_mat(mat_id))) {
		stop("cannot gen gstat for non existing mat ", mat_id)
	}
	downsample_n = scm_which_downsamp_n(scdb_mat(mat_id))
	gstat = scm_gene_stat(mat_id,
								 downsample_n = downsample_n)
	scdb_add_gstat(gstat_id, gstat)
}

#' Calculate basic statistics on a matrix
#'
#' @param mat The input matrix
#' @param niche_quantile A value between 0 and 1.
#'
#' @return Returns a dataset that contains statistic for all genes that passed the filtering stage.
#' Columns starting with ds contain UMI statistics after downsampling,
#' columns starting with n contain UMI statistics after normalizing UMIs so that the number of UMIs
#' per cell sums to 1.
#'
#' The columns are:
#'
#' \describe{
#'   \item{tot}{Total gene expression}
#'   \item{var}{Gene variance}
#'   \item{is_on_count}{Number of cells in which the gene is expressed}
#'   \item{sz_cor}{Correlation with cell size}
#'   \item{sz_cor_norm}{sz_cor after subtracting the trend}
#'   \item{niche_stat}{How many of the genen's umis are found in X\% of the most highly expressing cells. (regularized) }
#'   \item{niche_norm}{niche_stat after subtracting the niche_norm trend: median niche_norm value of genes with similar total expression}
#'   \item{n_mean}{Mean after normalization}
#'   \item{ds_top1}{Largest count, after downsampling}
#'   \item{ds_top2}{2nd largest count, after downsampling}
#'   \item{ds_top3}{3rd largest count, after downsampling}
#'   \item{ds_mean}{Mean on downsampled data}
#'   \item{ds_var}{Variance on downsampled data}
#'   \item{ds_log_varmean}{log2 of ds_var/ds_mean }
#'   \item{ds_vm_norm}{ds_log_varmean after subtracting the trend}
#'   \item{ds_is_on_count}{Number of cells in which the gene is expressed, after down sampling}
#'   \item{downsample_n}{Number of UMIs used for downsampling}
#'
#' }
# #'   \item{max_pk
# #' max_ratio
# #' var_meanpk
# #' ds_var_meanpk
#'
#'
#' @import zoo
#'
#' @export
#'
scm_gene_stat = function(mat_id,
                         niche_quantile = 0.2, #TODO: change to 0.1?
                         rseed=321, 
								 downsample_n = NULL,
								 K_std_n = 1000)
{
  oldseed = .set_seed(rseed)

#currently converting to non-sparse. Let's see if this need to be optimized
  scmat = scdb_mat(mat_id)
  mat =scmat@mat
  cat("Calculating gene statistics... ")

  if (niche_quantile >1 | niche_quantile < 0 ) {
    stop("niche_quantile must be between 0 and 1, got ", niche_quantile)
  }

  # filtering cells with too many or too few umis
	f_oversize = colSums(mat) > quantile(colSums(mat), 0.95) * 2 |
							colSums(mat) < 100
	
	# returns how many of the gene's umis are found in
	# X% of the most highly expressing cells. (regularized)
	quant_mean = function(x, k_reg) {
		n = length(x);
		up=sum(tail(sort(x), n=round(n*niche_quantile))); # sum umis in the top 20%
		return((k_reg+up)/(k_reg+sum(x)));
	}
	N = sum(!f_oversize)

	f_g = rowSums(mat) > 10

	downsample_n = scm_which_downsamp_n(scmat)

	message("will downsamp")
	mat_ds = scm_downsamp(mat, downsample_n)
	message("done downsamp")
	mat_ds = mat_ds[f_g,]
	mat_fg = mat[f_g,]

	message("will gen mat_n")
	mat_n = rescale_sparse_mat_cols(mat_fg, K_std_n/colSums(mat))
#	mat_n = t(t(mat_fg)*(K_std_n/colSums(mat)))
	message("done gen mat_n")

	n_ds = ncol(mat_ds)
	n = ncol(mat)
	stat_on_genes = function(i) {
		g_i = (1+(i-1)*quant):(min(nrow(mat_fg),i*quant))
		mat_gi = as.matrix(mat_fg[g_i,])
		mat_ds_gi = as.matrix(mat_ds[g_i,])
		gstat = data.frame(stringsAsFactors = F,
	   	name = rownames(mat_gi),
			tot = rowSums(mat_gi),
			var = round(matrixStats::rowVars(mat_gi),7),
			niche_stat =
		  		apply(mat_gi, 1, quant_mean, k_reg = 10),
			n_mean = round(matrixStats::rowMeans2(as.matrix(mat_n[g_i,])),7),
			ds_top1 = matrixStats::rowMaxs(mat_ds_gi) ,
			ds_top2 = matrixStats::rowOrderStats(mat_ds_gi, which = n_ds-1) ,
			ds_top3 = matrixStats::rowOrderStats(mat_ds_gi, which = n_ds-2) ,
			is_on_count = matrixStats::rowCounts(mat_ds_gi > 0),
			ds_is_on_count = matrixStats::rowCounts(mat_ds_gi > 0),
			ds_var = round(matrixStats::rowVars(mat_ds_gi),7),
			ds_mean = round(matrixStats::rowMeans2(mat_ds_gi),7))
		return(gstat)
	}
	max_bin = get_param("mc_cores")
	doMC::registerDoMC(max_bin)
	quant = ceiling(nrow(mat_fg)/max_bin)

	res <- plyr::alply(1:max_bin, 1, stat_on_genes, .parallel=TRUE)
	gene_stat = do.call(rbind, res)
	message("done computing basic gstat, will compute trends")

	if(ncol(mat) > 50000) {
		subs = sample(1:ncol(mat),50000)
		ctot = colSums(mat[,subs])
		gene_stat$sz_cor = round(apply(mat_fg[,subs], 1, function(x) { cor(x, ctot) }),3)

	} else {
		ctot = colSums(mat)
		gene_stat$sz_cor = round(apply(mat_fg, 1, function(x) { cor(x, ctot) }),3)
	}

	tot_ord = order(gene_stat$tot)
	cor_sz_ord = gene_stat$sz_cor[tot_ord]
	cmin = median(cor_sz_ord[1:101])
	cmax = median(cor_sz_ord[(length(cor_sz_ord)-101):length(cor_sz_ord)])
	sz_cor_trend = zoo::rollmedian(cor_sz_ord, 101, fill=c(cmin, NA, cmax))

	gene_stat$sz_cor_norm[tot_ord] =
		gene_stat$sz_cor[tot_ord] - sz_cor_trend

	# niche_stat =  how many of the genen's umis are found in
	# 20% of the most highly expressing cells. (regularized)
	niche_sz_ord = gene_stat$niche_stat[tot_ord]
	cmin = median(niche_sz_ord[1:101])
	cmax = median(niche_sz_ord[(length(niche_sz_ord)-101):length(niche_sz_ord)])
	niche_trend = zoo::rollmedian(niche_sz_ord, 101, fill=c(cmin, NA, cmax))

	gene_stat$niche_norm[tot_ord] =
		gene_stat$niche_stat[tot_ord] - niche_trend

	# var/mean
	ds_ord = order(gene_stat$ds_mean)

#	gene_stat$ds_log_varmean = ifelse(gene_stat$ds_mean>0.01, log2((0.001+gene_stat$ds_var)/(gene_stat$ds_mean+0.001)), 0)
	m_reg = 10/N
	gene_stat$ds_log_varmean = log2((m_reg+gene_stat$ds_var)/(m_reg+gene_stat$ds_mean))
	vm_sz_ord = gene_stat$ds_log_varmean[ds_ord]
	cmin = median(vm_sz_ord[1:101])
	cmax = median(vm_sz_ord[(length(vm_sz_ord)-101):length(vm_sz_ord)])
	vm_trend = zoo::rollmedian(vm_sz_ord, 101, fill=c(cmin, NA, cmax))

	gene_stat$ds_vm_norm[ds_ord] =
		gene_stat$ds_log_varmean[ds_ord] - vm_trend

	gene_stat$downsample_n = downsample_n

	# currently not used
	# gene_stat$max_ratio = gene_stat$max_pk/gene_stat$mean_pk
	# gene_stat$var_meanpk = gene_stat$var/gene_stat$mean_pk
	# gene_stat$ds_var_meanpk = gene_stat$ds_var/gene_stat$mean_pk

	rownames(gene_stat) = gene_stat$name

	gene_stat = gene_stat[,c("name","tot","var","is_on_count", "sz_cor","sz_cor_norm",
	                      "niche_stat",  "niche_norm", "n_mean",
	                      "ds_top1", "ds_top2", "ds_top3",
	                      "ds_mean", "ds_var", "ds_log_varmean","ds_vm_norm", "ds_is_on_count", "downsample_n")]

	.restore_seed(oldseed)
	cat("..done\n")

	return(gene_stat)
}


