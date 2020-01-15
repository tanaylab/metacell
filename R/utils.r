
# internal utility functions

# rowFunction is a function that works on rows, like rowMeans
# # much faster than tapply
.row_stats_by_factor = function (data, fact, rowFunction = rowMeans) {
	u = as.character(sort(unique(fact)))
	fact[is.na(fact)] = F
	n=length(u)
	centers = matrix(NA,dim(data)[1], n, dimnames = list(rownames(data), u))
	for (i in u) {
		if(sum(fact==i, na.rm=T)>1) {
			centers[,i] = rowFunction(data[,fact==i,drop=F])
		} else {
			centers[,i] = data[,fact==i]
		}
	} # much faster than tapply
	return(centers)
}

#' downsampl
#'
#' generate a umi matrix in which each cell (Column) is downsampled to a fixed number of molecules. This is a ehavy an inefficient proces at the moment, so we split it into chunks and run parallel using doMC
#'
#' @param umis the umi matrix
#' @param n the constant number of molecules per cell in the output matrix
#'
#' @export
#'
scm_downsamp = function(umis, n)
{
	old_seed = .set_seed(get_param("mc_rseed"))
	
	umis = umis[,colSums(umis)>= n]
	m = nrow(umis)
	.downsamp_one=function(v,n, replace = F){
		a = tabulate(sample(rep(1:length(v),times=v),replace=replace,size=n),nbins=m)
		return (a)
	}
	max_bin = get_param("mc_cores")
	doMC::registerDoMC(max_bin)

	max_bin = min(max_bin, ceiling(ncol(umis)/500))

	if(max_bin*10000 < ncol(umis)) {
		max_bin =  round(ncol(umis))/10000
	}
	cell_quant = ceiling(ncol(umis)/max_bin)
	seed = 19
	sub_dsamp = function(x) {
		set.seed(seed)
		i = 1+(x-1)*cell_quant
		j = min(x*cell_quant, ncol(umis))
	   ret = Matrix(apply(umis[,i:j], 2, .downsamp_one, n))
	   rownames(ret) = rownames(umis)
		return(as(ret,"dgCMatrix"))
	}
	res <- plyr::alply(1:max_bin, 1, sub_dsamp, .parallel=TRUE)
	umis_ds = do.call(cbind, res)

	.restore_seed(old_seed)

	return(umis_ds)
}

# save random number generator status, and set seed
.set_seed = function(rseed = 1) {
	if (exists(".Random.seed", .GlobalEnv)) {
		oldseed = .GlobalEnv$.Random.seed
	}
	else {
		oldseed = NULL
	}
	#message("will set seed")
	set.seed(seed=rseed)

	return(oldseed)
}

# restor random number generator status
.restore_seed = function(oldseed) {
	if (!is.null(oldseed)) {
		.GlobalEnv$.Random.seed = oldseed
	}
	else {
		rm(".Random.seed", envir = .GlobalEnv)
	}
}


# translate values vector to integers (for choosing colors by value from a color pallette
.vals_to_n = function(v, n = 1000, zlim = NULL) {

   if (is.null(zlim)) {
		zlim = range(v)
	}
	v = pmin(pmax(v, zlim[1]), zlim[2])
	v = v - min(v)
	1 + floor(v * ((n - 1) / max(v)))
}

# # from tlsrc/analysis/common/fread.r; one day this will be in a "common" package
# # returns rownames as the first column, named row.var
# ### @importFrom data.table fread
# fread <- function(...) data.table::fread(..., data.table=FALSE)
#
# fread_rownames <- function(..., row.var='rowname', set_rownames = F) {
# 	params <- list(...)
# 	header <- strsplit(readLines(params[[1]], n=1, warn=FALSE), '\t', fixed=TRUE)[[1]]
# 
# 	params$header = F
# 	params$skip = 1
# 	params$col.names = c(row.var, header)
# 	params$data.table = F
# 	
# 	mat = do.call(fread,params)
# 	
# 	if (set_rownames) {
# 		rownames(mat) = mat[,1]
# 		mat = mat[,-1]
# 	}
# 	return(mat)
# }

#' wrapping tgs functions to compute balanced graph from a matrix
#'
#' @param x matrix
#' @param knn K parameter
#' @param k_expand expansion parameter
#' @param k_alpha used to threshold blanaced scores (k_alpha*K > rank_in*rank_out0
#' @param k_beta used to threshold in-degree on k_beta*K before final thresholding of out-degree by K
#'
tgs_cor_graph = function(x, knn, k_expand, k_alpha, k_beta)
{
#	x_c = tgs_cor(x)
#	x_knn = tgs_knn(x_c, knn*k_expand)
	x_knn = tgs_cor_knn(x, y=x, knn=knn*k_expand)
	gr = tgs_graph(x_knn, knn, k_expand, k_beta)
	return(gr)
#	gr = tgs_cor_graph(x=feat, knn=K, k_expand=10, k_alpha=k_alpha, k_beta=k_beta)
}


#' Efficient version for scaling all columns in a sparse matrix
#'
#'
#' @param A first matrix
#' @param B second matrix - should be exactly same dimension as A
#'
#' @export

rescale_sparse_mat_cols = function(A, v_norm)
{
	if(!is(A, 'dgCMatrix')) {
		stop("cannot rescale sparse matrix for type that is not dgCMatrix")
	}
	A@x <- A@x * rep.int(v_norm, diff(A@p))
	return(A)
}

#' computing correlations between all rows in two matrices
#'
#' Efficient version for computing correlations between all rows
#'
#' @param A first matrix
#' @param B second matrix - should be exactly same dimension as A
#'
#' @export

allrow_cor = function(A,B) {
	cA <- A - rowMeans(A)
	cB <- B - rowMeans(B)
	sA <- sqrt(rowMeans(cA^2))
	sB <- sqrt(rowMeans(cB^2))

	rowMeans(cA * cB) / (sA * sB)
}
