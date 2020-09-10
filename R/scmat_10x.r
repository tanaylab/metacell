
#' Load a matrix from a 10x multi-batch dataset. The scdb version of mcell_read_multi_scmat_10x
#'
#' @param mat_nm - the name of the new matrix in scdb
#' @param dataset_table_fn - The index file of the 10x multi batch dataset. Fileds include: Batch.Set.ID, mat_fn, genes_fn, cells_fn. The first field must be specifieid for each row. The others may be NA, in which case the import process will search for base_dir/Batch.Set.ID/matrix.mtx, genes.tsv (or features.tsv), barcodes.tsv.
#' @param base_dir - directory where raw 10x files are located (dataset table points to subdirectories of this)
#' @param force - if true, will import from 10x files even when the matrix is present in the DB
#'
#' @export
#'
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
#' @param mat_nm  the name of the new matrix in scdb
#' @param base_dir a directory with data files - if this is specified, the matrix,gene and cells file names are determined by defulat
#' @param matrix_fn if base_dir is missing, this must define the matrix file name to be imported
#' @param genes_fn if base_dir is missing this must define the genes file name to be imported
#' @param cells_fn if base_dir is missing, this must define the cells file name to be imported
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
		if(!file.exists(genes_fn) & !grepl("^http", genes_fn)) {
			genes_fn = sprintf("%s/features.tsv", base_dir)
		}
		cells_fn = sprintf("%s/barcodes.tsv", base_dir)
	} else {
		if(is.null(matrix_fn)) {
			stop("MCERR - in importing a single 10x matrix, either specify a basedir or the individual matrix, gene and cells file")
		}
	}
	scdb_add_mat(mat_nm,
			scmat_read_scmat_10x(matrix_fn, genes_fn, cells_fn, dataset_id=mat_nm))
	return(TRUE)
}
# imports custom count matrix from 10x run. Matrix file must be in RDS format. Either sparse or regular matrix format are supported.
mcell_import_scmat_10x_custom = function(mat_nm,
		matrix_fn = NULL,
		force = FALSE)
{
	if(!scdb_is_valid()) {
		stop("MCERR - scdb is not initialized, cannot import")
	}
	if(!force & !is.null(scdb_mat(mat_nm))) {
		return(TRUE)
	}
	scdb_add_mat(mat_nm,
			scmat_read_scmat_10x_custom(matrix_fn, dataset_id=mat_nm))
	return(TRUE)
}

#' read multiple 10x umi matrices and merge them, based on a table defining the datasets. Field amp_batch_id from the table is added to the cell name to prevent cell names clashes.
#'
#' @param datasets_table_fn the index file of the 10x multi batch dataset. Fileds include: Batch.Set.ID, mat_fn, genes_fn, cells_fn. The first field must be specifieid for each row. The others may be NA, in which case the import process will search for base_dir/Batch.Set.ID/matrix.mtx, genes.tsv, barcodes.tsv.
#' @param base_dir where the input matrices reside.
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
			mat_fn = sprintf("%s/%s", base_dir, mat_fn)
		}
		genes_fn = dsets$genes_fn[i]
		if(is.na(genes_fn)) {
			genes_fn = sprintf("%s/%s/genes.tsv", base_dir, dnm)
		} else {
			genes_fn = sprintf("%s/%s", base_dir, genes_fn)
		}
		cells_fn = dsets$cells_fn[i]
		if(is.na(cells_fn)) {
			cells_fn = sprintf("%s/%s/barcodes.tsv", base_dir, dnm)
		} else {
			cells_fn = sprintf("%s/%s", base_dir, cells_fn)
		}
		amat = scmat_read_scmat_10x(matrix_fn = mat_fn,
										genes_fn = genes_fn,
										cells_fn = cells_fn,
										dataset_id = dnm)

		colnames(amat@mat) = paste(dnm, colnames(amat@mat), sep="_")
		amat@cells = colnames(amat@mat)
		rownames(amat@cell_metadata) = colnames(amat@mat)

		amat@cell_metadata$batch_set_id = dnm
		amat@cell_metadata$amp_batch_id = dnm
		amat@cell_metadata$seq_batch_id = dnm

		if(is.null(mat)) {
			message("will add")
			mat = amat
		} else {
			message("will merge")
			mat = scm_merge_mats(mat, amat)
		}
	}
	message("done reading")
	return(mat)
}


#' Read a matrix from the output of a 10x run. Batches can be stripped from the cell identifier if in BARCODE-LANE format.
#'
#' @param matrix_fn matrix input file name. Expecting a tab delimited file with 3 columns and a header: 'row', 'column', 'value'. row and column fields are the row and column matrix coordinates of the value.
#' @param genes_fn linking to the rownames.tab
#' @param cells_fn linking to the colnames.tab
#' @param paralog_policy - how to treat entries with the same gene name. Currently supporting "sum" for adding up the umis and "remove" for removing the gene. TBA: "split".
#'
#' @export
#'
#'
scmat_read_scmat_10x = function(matrix_fn,
		genes_fn,
		cells_fn,
		paralogs_policy = get_param("scm_10x_paralogs_policy"),
		use_batches = get_param("scm_10x_batches"),
		min_umis_n = get_param("scm_min_cell_umis"),
		dataset_id = "NS")
{
	remote_mode = F
	if(grepl("^http", matrix_fn)) {
		if(!url.exists(matrix_fn) |
			!url.exists(genes_fn) |
			!url.exists(cells_fn)) {
			stop("MC-ERR: missing 10x matrix urls, links: ", matrix_fn, " ", genes_fn, " ", cells_fn)
		} else {
			message("remote mode")
			remote_mode = T
		}
	} else if(!file.exists(matrix_fn) |
		!file.exists(genes_fn) |
		!file.exists(cells_fn)) {
		stop("MC-ERR: missing 10x matrix files, fns: ", matrix_fn, " ", genes_fn, " ", cells_fn)
	}
	# read sparse matrix with all genes and batches
	if(remote_mode) {
			download.file(matrix_fn, destfile="/tmp/matrix.mtx")
			download.file(genes_fn, destfile="/tmp/genes.tsv")
			download.file(cells_fn, destfile="/tmp/cells.tsv")
			matrix_fn = "/tmp/matrix.mtx"
			genes_fn = "/tmp/genes.tsv"
			cells_fn = "/tmp/cells.tsv"
	}
	umis = fread_mm(fname = matrix_fn, row.names = genes_fn, col.names = cells_fn)
	genes = as.data.frame(fread(genes_fn, header=F, stringsAsFactors = F))

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

	spike_regexp = get_param("scm_spike_regexp", 'metacell')
	if(!is.null(spike_regexp)) {
		ercc = grep(spike_regexp, rownames(umis))
		md$spike_count = colSums(umis[ercc,])
		no_ercc = grep(spike_regexp, rownames(umis), invert=T)
		umis = umis[no_ercc, ] # remove ERCC & first col
	} else {
		md$spike_count = 0
	}

	# filter small cells
	# f = colSums(umis) > min_umis_n
	# umis = umis[,f]
	# md = md[f,]
	# rownames(md) = colnames(umis)
	#
	# cat("MC: Filtered", sum(!f), "small cells\n")

	return(tgScMat(umis, stat_type="umi", cell_metadata = md))
}

#' Read a custom count matrix from the output of a 10x run. Batches can be stripped from the cell identifier if in BARCODE-LANE format.
#'
#' @param matrix_fn matrix input file name. Expecting a count matrix saved as .RDS file. Colnames represent barcodes, rownames represent custom gene annotation.
#'
#' @export
#'
#'
scmat_read_scmat_10x_custom = function(matrix_fn,
		#use_batches = get_param("scm_10x_batches"),
		#min_umis_n = get_param("scm_min_cell_umis"),
		dataset_id = "NS")
{
	if(!file.exists(matrix_fn)) 
	{
		stop("MC-ERR: missing matrix file: ", matrix_fn)
	}
	# read sparse matrix with all genes and batches
	umis = readRDS(matrix_fn)
	if (class(umis) == 'dgCMatrix')
	{
		umis = as.matrix(umis)
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
	# f = colSums(umis) > min_umis_n
	# umis = umis[,f]
	# md = md[f,]
	# rownames(md) = colnames(umis)
	#
	# cat("MC: Filtered", sum(!f), "small cells\n")

	return(tgScMat(umis, stat_type="umi", cell_metadata = md))
}
