#' Load a matrix from a simple dense table with genes in rows and cells in columns. First column is the gene name
#'
#' @param mat_nm - the name of the new matrix in scdb
#' @param fn  path of the matrix to read
#' @param meta_fn meta data table to read. If this is null, all cells will be onsidred as part of a single Amp.Batch and Dataset.ID
#' @param force - if true, will import from 10x files even when the matrix is present in the DB
#'
#' @export

mcell_import_scmat_tsv = function(mat_nm, fn, genes_fn = NULL, meta_fn = NULL, dset_nm = NULL, force = FALSE, paralogs_policy = "sum")
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
	umis = fread(fn)
	if(!is.null(genes_fn)) {
		all_ids = unlist(umis[,1])

		genes = as.data.frame(fread(genes_fn, header=F, stringsAsFactors = F))

		colnames(genes) = c('id', 'name')
		rownames(genes) = genes$id
		gname = genes[all_ids, 'name']
		umis = umis[!is.na(gname),]
		gname = gname[!is.na(gname)]

		unique_genes = names(which(table(gname) == 1))
		umis1 = umis[gname %in% unique_genes, -1]
		rownames(umis1) = gname[gname %in% unique_genes]
		if (paralogs_policy == 'remove') {
			umis = umis1
		}
		else if (paralogs_policy == 'sum') {
			non_unique_genes = setdiff(gname, unique_genes)
			umis2 = as.matrix(umis[gname %in% non_unique_genes, -1])

			message(sprintf("summing up total of %d paralog genes into %d unique genes", nrow(umis2), length(non_unique_genes)))

		   dup_genes = gname[gname %in% non_unique_genes]

			umis2s = apply(umis2, 2, function(x) {
					tapply(x, INDEX=dup_genes, FUN=sum) } )

			umis = rbind(umis1, umis2s)
			rownames(umis) = c(rownames(umis1), rownames(umis2s))
		}
	} else {
		rownames(umis) = umis[,1]
		umis = umis[,-1]
	}
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
				nrow=ncol(umis), ncol=4, byrow=T,
				dimnames=list(colnames(umis),
							c('type', 'batch_set_id', 'amp_batch_id', 'seq_batch_id')))
		md = as.data.frame(md)
	}
	sparse_m = Matrix(as.matrix(umis))
	rownames(sparse_m) = rownames(umis)
	scmat = scm_new_matrix(sparse_m, md)
	scdb_add_mat(mat_nm, scmat)
	return(TRUE)
}
