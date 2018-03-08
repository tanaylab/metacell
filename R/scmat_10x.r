
#' Load a matrix from a 10x multi-batch dataset. The scdb version of mcell_read_multi_scmat_10x
#'
#' @param mat_nm - the name of the new matrix in scdb
#' @param dataset_table_fn - path of the key table.
#' @param base_dir - directory where raw 10x files are located (dataset table points to subdirectories of this)
#' @param force - if true, will import from 10x files even when the matrix is present in the DB
#'
#' @export

mcell_import_multi_scmat_10x = function(mat_nm, 
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
	scdb_add_mat(mat_nm, mcell_read_multi_scmat_10x(dataset_table_fn, base_dir))
	return(TRUE)
}

#' Load a matrix from a 10x dataset. The scdb version of mcell_read_multi_scmat_10x
#'
#' @param mat_nm - the name of the new matrix in scdb
#' @param dataset_table_fn - path of the key table
#' @param force - if true, will import from 10x files even when the matrix is present in the DB
#'
#' @export

mcell_import_scmat_10x = function(mat_nm, 
		base_dir = NULL, 
		matrix_fn = NULL, 
		genes_fn = NULL,
		cells_fn = NULL,
		force = FALSE)
{
	if(!scdb_is_valid()) {
		stop("MCERR - scdb is not initialized, cannot import")
	}
	if(!force & !is.null(scdb_mat(mat_nm))) {
		return(TRUE)
	}
	if(!is.null(base_dir)) {
		matrix_fn = sprintf("%s/matrix.mtx", base_dir)
		genes_fn = sprintf("%s/genes.tsv", base_dir)
		cells_fn = sprintf("%s/barcodes.tsv", base_dir)
	}
	scdb_add_mat(mat_nm, 
			mcell_read_scmat_10x(matrix_fn, genes_fn, cells_fn, dataset_id=mat_nm))
	return(TRUE)
}

#' read multiple 10x umi matrices and merge them, based on a table defining the datasets.
#' 
#' @param datasets_table_fn the index file of the 10x multi batch dataset. Fileds include: Batch.Set.ID, mat_fn, genes_fn, cells_fn. The first field must be specifieid for each row. The others may be NA, in which case the import process will search for base_dir/Batch.Set.ID/matrix.mtx, genes.tsv, barcodes.tsv. 
#' @param mcell_min_cell_umis (ENVIRON) minimum number of umi's for retaining the cell in the matrix (default: 200)
#'
#' @export

mcell_read_multi_scmat_10x = function(datasets_table_fn, base_dir)
{
	if(!file.exists(datasets_table_fn)) {
		stop("MC-ERR: 10x multi_batch index file is missing, fn = ", datasets_table_fn)
	}
	dsets = fread(datasets_table_fn)
	miss_f = setdiff(c("Batch.Set.ID", "mat_fn", "genes_fn", "cells_fn" ),
					colnames(dsets))
	if(length(miss_f)>0) {
		stop("MC-ERR: missing fields in 10x dataset table fn : ", miss_f)
	}

	mat = NULL

	for(i in 1:nrow(dsets)) {
		dnm = dsets$Batch.Set.ID[i]
		message("will read ", dnm)
		mat_fn = dsets$mat_fn[i]
		if(is.na(mat_fn)) {
			mat_fn = sprintf("%s/%s/matrix.mtx", base_dir, dnm)
		} else {
			mat_fn = sptrinf("%s/%s", base_dir, mat_fn)
		}
		genes_fn = dsets$genes_fn[i]
		if(is.na(genes_fn)) {
			genes_fn = sprintf("%s/%s/genes.tsv", base_dir, dnm)
		} else {
			genes_fn = sptrinf("%s/%s", base_dir, genes_fn)
		}
		cells_fn = dsets$cells_fn[i]
		if(is.na(cells_fn)) {
			cells_fn = sprintf("%s/%s/barcodes.tsv", base_dir, dnm)
		} else {
			cells_fn = sptrinf("%s/%s", base_dir, cells_fn)
		}
		amat = mcell_read_scmat_10x(matrix_fn = mat_fn,
										genes_fn = genes_fn,
										cells_fn = cells_fn,
										dataset_id = dnm)
		if(is.null(mat)) {
			mat = amat
			mat@cell_metadata$batch_set_id = dnm
			mat@cell_metadata$amp_batch_id = dnm
			mat@cell_metadata$seq_batch_id = dnm
		} else {
			mat = scm_merge_mats(mat, amat, dnm)
		}
	}
	return(mat)
}


#' Read a matrix from the output of a 10x run. Batches can be stripped from the cell identifier if in BARCODE-LANE format.
#'
#' @param matrix_fn matrix input file name. Expecting a tab delimited file with 3 columns and a header: 'row', 'column', 'value'. row and column fields are the row and column matrix coordinates of the value. 
#' @param genes_fn linking to the rownames.tab
#'	@param cells_fn linking to the colnames.tab 
#' @param min_umis_n minimum number of umi's for retaining the cell in the matrix (default: 200)
#' @param paralog_policy - how to treat entries with the same gene name. Currently supporting "sum" for adding up the umis and "remove" for removing the gene. TBA: "split".
#'
#' @export
#'
#' @importFrom data.table fread
#'
mcell_read_scmat_10x = function(matrix_fn, 
		genes_fn,
		cells_fn,
		paralogs_policy = get_param("scm_10x_paralogs_policy"),
		use_batches = get_param("scm_10x_batches"),
		min_umis_n = get_param("scm_min_cell_umis"),
		dataset_id = "NS")
{
	if(!file.exists(matrix_fn) |
		!file.exists(genes_fn) |
		!file.exists(cells_fn)) {
		stop("MC-ERR: missing 10x matrix files, fns: ", matrix_fn, " ", genes_fn, " ", cells_fn)
	}
	# read sparse matrix with all genes and batches
	umis = fread_mm(fname = matrix_fn, row.names = genes_fn, col.names = cells_fn)
	
	genes = read.table(genes_fn, header=F, stringsAsFactors = F)
	colnames(genes) = c('id', 'name')
	rownames(genes) = genes$id
	
	unique_genes = names(which(table(genes$name) == 1))
	umis1 = umis[genes$name %in% unique_genes, ]
	rownames(umis1) = genes[rownames(umis1), 'name']
	if (paralogs_policy == 'remove') {
		umis = umis1
	}
	else if (paralogs_policy == 'sum') {
		non_unique_genes = setdiff(genes$name, unique_genes)
		umis2 = as.matrix(umis[genes$name %in% non_unique_genes, ])
		
		message(sprintf("summing up total of %d paralog genes into %d unique genes", nrow(umis2), length(non_unique_genes)))
		
		umis2s = apply(umis2, 2, function(x) { tapply(x, INDEX=genes[rownames(umis2), 'name'], FUN=sum) } )
		
		umis = rbind(umis1, umis2s)
	}
	else {
		stop(sprintf("MC-ERR: Loading 10x data, unknown paralogs policy (%s), supprting remove/sum", paralogs_policy))
	}
	
	cell_lane = gsub(".*-", "", colnames(umis))
	if(max(length(cell_lane)) == 0) {
		cell_batch = dataset_id
	} else {
		cell_lane[cell_lane == ""] = 0
		cell_batch=paste(sprintf("%s_L", dataset_id), cell_lane, sep="")
	}
	
	md = data.frame(row.names=colnames(umis), 
						 type='10x', 
						 batch_set_id=cell_batch,
						 amp_batch_id=cell_batch,
						 seq_batch_id=cell_batch)

	spike_regexp = get_param("scm_spike_regexp")
	if(!is.null(spike_regexp)) {	
		ercc = grep(spike_regexp, rownames(umis))
		md$spike_count = colSums(umis[ercc,])
		no_ercc = grep(spike_regexp, rownames(umis), invert=T)
		umis = umis[no_ercc, ] # remove ERCC & first col
	} else {
		md$spike_count = 0
	}
	
	# filter small cells
	f = colSums(umis) > min_umis_n
	umis = umis[,f]
	md = md[f,]
	rownames(md) = colnames(umis)

	cat("MC: Filtered", sum(!f), "small cells\n")

	return(tgScMat(umis, stat_type="umi", cell_metadata = md))			
}

