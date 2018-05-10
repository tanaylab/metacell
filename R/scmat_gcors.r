
#' Find and report naive gene gene correlation over umi matrix
#'
#' The generated table should be the basis for manual curation usually, and can be used to generate supervised gene sets, or a basis for them
#'
#' @param mat_id matrix to use
#' @param gene_anchors genes that will be used as anchors for computing correlation against all other genes
#' @param gene_anti genes that will be used as negative anchors - high correlation to these will be considered negatively
#' @param cor_thresh correlation threshold (on downsampled data) - the table will include genes with maximal correlation higher than the threshold, AND, maximal correlation for negative anchor lower than for positive nachors
#' @param tab_fn file name of tabular report
#' @export
#'

mcell_mat_rpt_cor_anchors = function(mat_id, gene_anchors, gene_anti = c(), cor_thresh, tab_fn, sz_cor_thresh = NA, downsample_n = NA, cells = NULL)
{
	mat = scdb_mat(mat_id)
	if(is.null(mat)) {
		stop("missing mat ", mat_id)
	}

	if(!is.null(cells)) {
		umis = mat@mat[,cells]
	} else {
		umis = mat@mat
	}

	if(is.na(downsample_n)) {
		downsample_n = quantile(colSums(umis), 0.05)
	}
	
	mat_ds = scm_downsamp(umis, downsample_n)
	mat_ds = mat_ds[rowSums(mat_ds)>10,]
	csize = colSums(mat@mat[,colnames(mat_ds)])
	gcors = data.frame(sz_cor = apply(mat_ds, 1, cor, csize))
	for(g in gene_anchors) {
		gcors[,g] = apply(mat_ds, 1, cor, mat_ds[g,])
	}
	for(g in gene_anti) {
		gcors[,g] = apply(mat_ds, 1, cor, mat_ds[g,])
	}
	N1 = length(gene_anchors) + 1
	N2 = N1 + length(gene_anti)

	if(length(gene_anchors) == 1) {
		gcors$max = gcors[,2]
	} else {
		gcors$max = apply(gcors[,2:N1], 1, max)
	}
	if(length(gene_anti) == 0) {
		gcors$neg_max = 0
	} else if(length(gene_anti) == 1) {
		gcors$neg_max = gcors[,N2]
	} else {
		gcors$neg_max = apply(gcors[,(N1+1):N2], 1, max)
	}
	f = !is.na(gcors$max) & gcors$max > cor_thresh & gcors$max > gcors$neg_max
	if(!is.na(sz_cor_thresh)) {
		f = f | (gcors$sz_cor > sz_cor_thresh)
	}
	write.table(gcors[f ,], tab_fn, sep="\t", quote=F)
}
