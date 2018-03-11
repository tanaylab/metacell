#' Load a matrix from a MARS-seq multi-batch dataset. The scdb version of mcell_read_multi_scmat_mars
#'
#' @param mat_nm - the name of the new matrix in scdb
#' @param dataset_table_fn - path of the key table.
#' @param base_dir - directory where raw mars files are located (dataset table points to subdirectories of this)
#' @param force - if true, will import from MARS files even when the matrix is present in the DB
#'
#' @export

mcell_import_multi_mars = function(mat_nm,
							dataset_table_fn,
							base_dir,
							force=FALSE)
{
	if(!scdb_is_valid()) {
		stop("MCERR - scdb is not initialized, cannot import")
	}
	if(!force & !is.null(scdb_mat(mat_nm))) {
		return(TRUE)
	}
	scdb_add_mat(mat_nm,
			 mcell_read_multi_scmat_mars(dataset_table_fn, base_dir))
	return(TRUE)
}

#' read multiple MARS umi matrices and merge them, based on a table defining the datasets.
#'
#' @param datasets_table_fn the index file of the MARS multi batch dataset. This is a tab delimited text file, with an arbitrary number of columns and a header line. The three mandatory fields are Amp.Batch.ID and Seq.Batch.ID, Batch.Set.ID. The first specify the ID of the batch defined by the row, and also the file name (without the .txt suffix) of the respective umi table in the base_dir provided. The second defines and ID of the sequencing batch (may be relevant for further noise cleanups beyond those done in the low-level pipeline). The third id group different batches into sets for downstream analysis (e.g. QC and more). This is specifically geared toward the output of the std MARS pipeline circa 2014...(bummer).
#' @param base_dir defines the umitab directory
#'
#' @export

mcell_read_multi_scmat_mars = function(datasets_table_fn, base_dir)
{
	if(!file.exists(datasets_table_fn)) {
		stop("MC-ERR: MARS multi_batch index file is missing, fn = ", datasets_table_fn)
	}
	dsets = fread(datasets_table_fn, sep="\t")
	mandatory = c("Amp.Batch.ID", "Seq.Batch.ID", "Batch.Set.ID")
	miss_f = setdiff(mandatory, colnames(dsets))
	if(length(miss_f)>0) {
		stop("MC-ERR: missing fields in MARS dataset table fn : ", miss_f)
	}

	mat = NULL

	for(i in 1:nrow(dsets)) {
		amp_batch = dsets$Amp.Batch.ID[i]
		message("will read ", amp_batch)

		umis = fread_rownames(sprintf("%s/%s.txt", base_dir, amp_batch), sep="\t", set_rownames=T)

		md = as.data.frame(matrix(c('MARS', unlist(dsets[i,])), nrow=ncol(umis), ncol=ncol(dsets)+1, byrow=T, dimnames=list(colnames(umis), c('type', colnames(dsets))))) %>%
		  rename(batch_set_id=Batch.Set.ID, amp_batch_id=Amp.Batch.ID, seq_batch_id=Seq.Batch.ID)

		spike_regexp = get_param("scm_spike_regexp")
		if(!is.null(spike_regexp)) {
		  ercc = grep(spike_regexp, rownames(umis))
		  md$spike_count = colSums(umis[ercc,])
		  no_ercc = grep(spike_regexp, rownames(umis), invert=T)
		  umis = umis[no_ercc, ] # remove ERCC & first col
		} else {
		  md$spike_count = 0
		}

		amat = tgScMat(as.matrix(umis), stat_type = "umi", cell_metadata = md)

		if(is.null(mat)) {
			mat = amat
		} else {
			mat = scm_merge_mats(mat, amat)
		}

	}
	return(mat)
}

