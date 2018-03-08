
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
	umis = umis[,colSums(umis)>= n]
	m = nrow(umis)
	.downsamp_one=function(v,n, replace = F){
		a = tabulate(sample(rep(1:length(v),times=v),replace=replace,size=n),nbins=m)
		return (a)
	}
	max_bin = get_param("mc_cores")
	doMC::registerDoMC(max_bin)

	if(max_bin*10000 < ncol(umis)) {
		max_bin =  round(ncol(umis))/10000
	}
	cell_quant = ceiling(ncol(umis)/max_bin)
	sub_dsamp = function(x) {
		i = 1+(x-1)*cell_quant 
		j = min(x*cell_quant, ncol(umis))  
	   ret = Matrix(apply(umis[,i:j], 2, .downsamp_one, n))
	   rownames(ret) = rownames(umis)
		return(as(ret,"dgCMatrix"))
	}
	res <- plyr::alply(1:max_bin, 1, sub_dsamp, .parallel=TRUE)
	umis_ds = do.call(cbind, res)
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
	set.seed(rseed)

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

# from tlsrc/analysis/common/fread.r; one day this will be in a "common" package
# returns rownames as the first column, named row.var
### @importFrom data.table fread
fread <- function(...) data.table::fread(..., data.table=FALSE)

fread_rownames <- function(..., row.var='rowname', set_rownames = F) {
	params <- list(...)
	header <- strsplit(readLines(params[[1]], n=1, warn=FALSE), '\t', fixed=TRUE)[[1]]

	params$header = F; params$skip = 1; params$col.names = c(row.var, header)

	mat = do.call(fread,params)
	if (set_rownames) {
		rownames(mat) = mat[,1]
		mat = mat[,-1]
	}
	return(mat)
}


