#' Load a matrix from a simple dense table with genes in rows and cells in columns. First column is the gene name
#'
#' @param mat_nm - the name of the new matrix in scdb
#' @param fn  path of the matrix to read
#' @param meta_fn meta data table to read. If this is null, all cells will be onsidred as part of a single Amp.Batch and Dataset.ID
#' @param force - if true, will import from 10x files even when the matrix is present in the DB
#'
#' @export

mcell_import_scmat_tsv = function(mat_nm, fn, meta_fn = NULL, dset_nm = NULL, force = FALSE)
{
	if(!scdb_is_valid()) {
		stop("MCERR - scdb is not initialized, cannot import")
	}
	if(!force & !is.null(scdb_mat(mat_nm))) {
		return(TRUE)
	}
	if(!file.exists(fn)) {
		stop("MCERR - umi matrix file ", fn, " does not exist")
	}
	mat = fread(fn)
	rownames(mat) = mat[,1]
	mat = mat[,-1]
	if(!is.null(meta_fn)) {
		md = fread(meta_fn, sep="\t")
		mandatory = c("Amp.Batch.ID", "Seq.Batch.ID", "Batch.Set.ID")
		miss_f = setdiff(mandatory, colnames(md))
		if(length(miss_f)>0) {
			stop("MC-ERR: missing fields in matrix metadata : ", miss_f)
		}
	} else {
		if(is.null(dset_nm)) {
			stop("either provide a meta data file, or supply a dataset name (parameter dset_nm) and we'll fake it for you...")
		}
		dsets = data.frame(batch_set_id = dset_nm,
							 	 amp_batch_id = dset_nm,
							    seq_batch_id = dset_nm)
		md = matrix(c('import', dset_nm, dset_nm, dset_nm), 
				nrow=ncol(mat), ncol=4, byrow=T, 
				dimnames=list(colnames(mat), 
							c('type', 'batch_set_id', 'amp_batch_id', 'seq_batch_id')))
		md = as.data.frame(md)
	}
	scmat = scm_new_matrix(Matrix(as.matrix(mat)), md)
	scdb_add_mat(mat_nm, scmat)
	return(TRUE)
}
