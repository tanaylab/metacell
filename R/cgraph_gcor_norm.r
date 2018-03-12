#' Computing gene-gene correlation normaized over a similarity graph
#'
#' This is generating a randomized umi matrix by substituting each umi call (per cell and gene) with a umi call from the same gene on a cell that is connected to the original cell in a given similarity graph. If the similarity graph is based on some specific features that express a specific process (e.g. cell cycle), this can help in normalizing this effect and retaining correlation that are independent of it.
#'
#' @param cgraph_id id of the cell graph to be used for normalization
#' @param mat_id the matrix with umis to be analyzed
#' @param K the maxium number of edges to consider per cell in normalization. This can be smaller than the K used in cgraph_id, to provide tighter control
#' @param min_gtot minimal number of umis to consider when computing gene-gene correlaitons
#' @return a list with two entries. gcor - is the gene correlation matrix of the downsampled matrix. r_gcor - is the correlation of the randomized umi matrix
#'
#' @export

mcell_cgraph_norm_gcor = function(cgraph_id, mat_id, K=-1, min_gtot=1000)
{
	cgraph = scdb_cgraph(cgraph_id)
	if(is.null(cgraph)) {
		stop("missing cgraph for gene cor normalization, id = ", cgraph_id)
	}
	mat = scdb_mat(mat_id)
	if(is.null(mat)) {
		stop("missing mat for gene cor normalization, id = ", mat_id)
	}

	dsamp_n = scm_which_downsamp_n(mat)

	mat_ds = scm_downsamp(mat@mat, dsamp_n)
	c_nms = colnames(mat_ds)

	adjs = split(cgraph@edges$mc2, cgraph@edges$mc1)
	adjs = adjs[c_nms]
	if(K==-1) {
	  K = max(unlist(lapply(adjs, length)))
	}
	m_adjs = t(sapply(adjs, function(x) {
						x1 = intersect(x, c_nms);
						x[1+(0:(K-1))%%length(x)]
						}))
	C = nrow(m_adjs)
	f_gene = rowSums(mat_ds)>min_gtot
	r_mat = t(apply(mat_ds[f_gene,], 1, function(x) { 
				shuf = m_adjs[(1:C)+C*(sample(1:K,s=C, r=T)-1)]
				x[shuf]
				}))

	gcor = cor(as.matrix(t(mat_ds[f_gene,])), m="spearman")
	r_mat[is.na(r_mat)] = 0
	r_gcor = cor(t(r_mat), m="spearman")
	return(list(gcor=gcor, r_gcor=r_gcor))
}
