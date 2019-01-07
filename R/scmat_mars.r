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
							patch_cell_name=F,
							force=FALSE)
{
	if(!scdb_is_valid()) {
		stop("MCERR - scdb is not initialized, cannot import")
	}
	if(!force & !is.null(scdb_mat(mat_nm))) {
		return(TRUE)
	}
	scdb_add_mat(mat_nm,
			 mcell_read_multi_scmat_mars(dataset_table_fn, base_dir, patch_cell_name=patch_cell_name))
	return(TRUE)
}

#' read multiple MARS umi matrices and merge them, based on a table defining the datasets.
#'
#' @param datasets_table_fn the index file of the MARS multi batch dataset. This is a tab delimited text file, with an arbitrary number of columns and a header line. The three mandatory fields are Amp.Batch.ID and Seq.Batch.ID, Batch.Set.ID. The first specify the ID of the batch defined by the row, and also the file name (without the .txt suffix) of the respective umi table in the base_dir provided. The second defines and ID of the sequencing batch (may be relevant for further noise cleanups beyond those done in the low-level pipeline). The third id group different batches into sets for downstream analysis (e.g. QC and more). This is specifically geared toward the output of the std MARS pipeline circa 2014...(bummer).
#' @param base_dir defines the umitab directory
#'
#' @export

mcell_read_multi_scmat_mars = function(datasets_table_fn, base_dir, patch_cell_name=F)
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

		umis = fread_rownames(sprintf("%s/%s.txt", base_dir, amp_batch), sep="\t", row.var=NULL)
		if(patch_cell_name) {
			cn = paste(amp_batch, colnames(umis), sep=".")
			colnames(umis) = cn
		}

		md = as.data.frame(matrix(c('MARS', unlist(dsets[i,])), nrow=ncol(umis), ncol=ncol(dsets)+1, byrow=T, dimnames=list(colnames(umis), c('type', colnames(dsets)))))
		md = rename(md, batch_set_id=Batch.Set.ID, amp_batch_id=Amp.Batch.ID, seq_batch_id=Seq.Batch.ID)

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

#' Creates a new matrix object from a given one by adding FACS index sorting data to matrix metadata table
#'
#'#' @param mat_nm - the name of the new matrix in scdb
#' @param dataset_table_fn - path of the key table.
#' @param base_dir - directory where raw mars files are located (dataset table points to subdirectories of this)

#' @param new_mat_id - id of the new matrix in scdb
#' @param mat_id - id of the original matrix to add data to
#' @param base_dir - directory where raw tab-delimited index-sort data files are located. Expecting file name to be <amp_batch_id>.txt with well names as row names and antibodies as column names
#' @param amp_batch_ids - if not null, upload data from amp_batch_id specified here.
#' @param force overwrite output matrix if it exist
#'
#' @export
#'

mcell_add_mars_facs_data = function(new_mat_id, mat_id,
                                   base_dir, amp_batch_ids = NULL, force=T)
{
  if(!scdb_is_valid()) {
    stop("MCERR - scdb is not initialized, cannot add data")
  }
  if(!force & !is.null(scdb_mat(mat_id))) {
    return(TRUE)
  }
  scdb_add_mat(new_mat_id,
               scm_add_mars_facs_data(mat_id, base_dir, amp_batch_ids))
  return(TRUE)
}

#' Adds FACS index sorting data to matrix metadata table
#'
#' @param mat_id - id of the original matrix to add data to
#' @param base_dir - directory where raw tab-delimited index-sort data files are located. Expecting file name to be <amp_batch_id>.txt with well names as row names and antibodies as column names
#' @param amp_batch_ids - if not null, upload data from amp_batch_id specified here.
#'
#' @return new matrix with updated metadata table
#'
scm_add_mars_facs_data = function(mat_id, base_dir, amp_batch_ids = NULL)
{
	scmat = scdb_mat(mat_id) 
	if (is.null(scmat)) {
		stop(sprintf("MC-ERR: trying to add FACS data to an undefined mat object %s", mat_id))
	}
  wells2cells_fn = get_param("scm_mars_wells2cells_fn")
  if (is.null(wells2cells_fn) | !file.exists(wells2cells_fn)) {
    stop("MC-ERR: Missing wells2cells file when adding MARS FACS data")
  }
  w2c_batch = get_param("scm_mars_wells2cells_batch_field")
  w2c_cell = get_param("scm_mars_wells2cells_cell_field")
  w2c_well = get_param("scm_mars_wells2cells_well_field")

  all_batches = unique(scmat@cell_metadata[, "amp_batch_id"])

  if (is.null(amp_batch_ids)) {
    amp_batch_ids = all_batches
  } else {
    dif = setdiff(amp_batch_ids, all_batches)
    if (length(dif > 0)) {
      stop ("MC-ERR: FACS data tables supplied for batches not loaded: ", paste0(dif, collapse=", "))
    }
  }

  wells_dict = read.table(wells2cells_fn, header=T, sep="\t", stringsAsFactors = F)
  rownames(wells_dict) = wells_dict[, w2c_cell]

  missing_cells = setdiff(rownames(scmat@cell_metadata), rownames(wells_dict))
  if (length(missing_cells) > 0) {
    stop(sprintf("MC-ERR: some cells in mat are not found in given wells_cells file (%s), first 10: %s", wells2cells_fn, paste0(missing_cells[1:10], collapse = ", ")))
  }
  wells_dict = wells_dict[ rownames(scmat@cell_metadata), c(w2c_batch, w2c_cell, w2c_well)]

  idx_tab = NULL
  for (batch in amp_batch_ids) {
    ifn = sprintf("%s/%s.txt", base_dir, batch)
    if (file.exists(ifn)) {
      facs_idx = read.table(ifn,
                            header = T,
                            sep = "\t",
                            stringsAsFactors = F)
      facs_idx[, w2c_batch] = batch

      if (is.null(idx_tab)) {
        idx_tab = facs_idx
      }
      else {
        idx_tab[, setdiff(colnames(facs_idx), colnames(idx_tab))] = NA
        facs_idx[, setdiff(colnames(idx_tab), colnames(facs_idx))] = NA
        idx_tab = rbind(idx_tab, facs_idx[, colnames(idx_tab)])
      }
    }
  }
  message("loaded facs idxs from batches, merging with wells_cells...")
  m = merge(wells_dict, idx_tab, by.x=c(w2c_batch, w2c_well), by.y=c(w2c_batch, "Well"), all.x=T)
  rownames(m) = m[, w2c_cell]
  colnames(m) = paste0(colnames(m), "_Ab")
  scmat@cell_metadata = cbind(scmat@cell_metadata, m[ rownames(scmat@cell_metadata), ])

  scmat
}
