# this is important for running this package from scripts
#' @import methods
NULL

#' Single cell RNA-seq matrix
#'
#' a single cell RNA matrix interface. Does not do much beyond adding
#' specifically gene names, cell names and labels describing the type of data.
#' Methods for basic stat, gene selection, clustering are implemented separately
#'
#' @slot mat A sparse matrix containing the expression data
#' @slot genes list of all genes. By default the matrix row names are used
#' @slot cells list of all cells. By default the matrix col names are used
#' @slot stat_type Type of statistic. Currently only "umi" is supported
#' @slot cell_metadata dataframe with metadata on the cells. Rownames correspond to the cell names. By default a table consisting of only a "Batch" field is added with a trivial batch number. If the cell_metadata does not include a Batch field, it will be added with a trivial value. Additional standard fields spike_umis - the number of (filtered_) spike umis for the cell
#'
#' @examples
#' # we first initialize a tgScMat object
#' #TO ADD
#' @importClassesFrom Matrix Matrix dgCMatrix dgeMatrix
#' @import Matrix
#'
#' @export tgScMat
#' @exportClass tgScMat
tgScMat <- setClass(
   "tgScMat",
	slots = c(
	  mat = "dgCMatrix",
	  genes = "vector",
	  cells = "vector",
	  ncells = "numeric",
	  ngenes = "numeric",
	  stat_type = "character",
	  ignore_genes = "vector",
	  ignore_cells = "vector",
	  ignore_gmat = "dgCMatrix",
	  ignore_cmat = "dgCMatrix",
	  ignore_gcmat = "dgCMatrix",
	  cell_metadata = "data.frame")
)

setMethod(
  "initialize",
  signature = "tgScMat",
  definition =
    function(.Object, mat=NULL, stat_type = "umi", cell_metadata = NULL, ...) {
      .Object@stat_type=stat_type
      if(!is.null(mat)) {
      	if(!is.null(rownames(mat))) {
        		.Object@genes=rownames(mat)
      	} else {
        		.Object@genes=seq(1,nrow(mat),1)
      	}
      	if(!is.null(colnames(mat))) {
        		.Object@cells=colnames(mat)
      	} else {
      	  .Object@cells=seq(1,ncol(mat),1)
      	}
			.Object@ngenes = nrow(mat)
      	.Object@ncells = ncol(mat)

      	.Object@mat = as(mat, "dgCMatrix")

      	# cell metadata
      	.Object@cell_metadata = .scm_init_md(.Object, cell_metadata)
      }

      return(.Object)
    }
)

# Checks metadata
.scm_init_md = function(.Object, cell_metadata)
{
  if(is.null(cell_metadata)) {
		#by defulat the mat data indicate one single Batch
		md = data.frame(amp_batch_id = rep(1, nrow(.Object@mat)),
							 row.names=.Object@cells)
		return(md)
  }
  md_nms = rownames(cell_metadata)
  if(length(md_nms) != length(.Object@cells)
	  | length(md_nms) != length(intersect(.Object@cells, md_nms))) {
		stop("Metadata cells names are incompatible with matrix cells - aborting")
  }
  if(!"amp_batch_id" %in% colnames(cell_metadata)) {
		cell_metadata$amp_batch_id = 1
  }
  cell_metadata = cell_metadata[.Object@cells,]	#making sure it is reordered

  return(cell_metadata)
}

setValidity("tgScMat", function(object) {
  # TBA
  return(TRUE)
})



# PUBLIC CONSTRUCTORS

#' Constract a tgScMat
#'
#' Constract a tgScMat gene expression matrix object
#'
#'
#' @param mat A matrix containing the expression data (will be stored as a sparse matrix)
#' @param cell_metadata dataframe with metadata on the cells. Rownames correspond to the cell names.
#' @param stat_type Type of statistic. ("umi", "rpkm", "tpm" etc.). Default: "umi"

#'
#' @return A tgScMat gene expression matrix object
#'
#' @export
scm_new_matrix = function(mat,
									cell_metadata,
									stat_type = "umi")
{
	return(
		tgScMat(
			mat,
			stat_type = stat_type,
			cell_metadata = cell_metadata
		)
	)
}

#' ncell_merge_mats: Merge two matrix using their ids in scdb. See scm_merge_mats for details on batch management and policies on missing genes
#'
#' @param id1 - first matrix to merge
#' @param id2 - second matrix to merge
#' @param new_id - id of matrix to generate

mcell_merge_mats = function(id1, id2, new_id)
{
	if(!is.character(id1) | !is.character(id2) | !is.character(new_id)) {
		stop("mcell_merge_mats must be called ids into scdb, not other objects")
	}
	if(is.null(scdb_mat(id1)) | is.null(scdb_mat(id2))) {
		stop("mcell_merge_mats called with missing mat ids, ", id1, " ", id2)
	}
	scdb_add_mat(new_id,
					 scm_merge_mats(scdb_mat(id1), scdb_mat(id2)))
}

#' scm_merge_mats: Merge two single cell matrix object. Return the merged matrix, with merged meta data and issues an error if there are overlapping cell names between the two matrices. In case genes sets differs between the matrices, the union is used, with zeros (not NAs!) filling up the missing genes in the respective matrix.
#'
#' @param scmat1  tgScMat object.
#' @param scmat2  tgScMat object.
#'
#' @export
#'
scm_merge_mats = function(scmat1, scmat2)
{
	if(class(scmat1)[1] != "tgScMat") {
		stop("invalid scmat1 in scm_merge_mat - if you want to merge from the scdb, call mcell_merge_mats")
	}
	if(class(scmat2)[1] != "tgScMat") {
		stop("invalid scmat2 in scm_merge_mat - if you want to merge from the scdb, call mcell_merge_mats")
	}

	scmat1 = scm_ignore_cells(scmat1, NULL)
	scmat2 = scm_ignore_cells(scmat2, NULL)
	md1 = scmat1@cell_metadata
	md2 = scmat2@cell_metadata
	md_f1 = colnames(md1)
	md_f2 = colnames(md2)
	if(length(setdiff(md_f2,md_f1))>0) {
		md1[,setdiff(md_f2,md_f1)]=NA
	}
	if(length(setdiff(md_f1,md_f2))>0) {
		md2[,setdiff(md_f1,md_f2)]=NA
	}
	scmat1@cell_metadata = rbind(md1, md2)
#	scmat1@cell_metadata = bind_rows(scmat1@cell_metadata, scmat2@cell_metadata)

	m1 = scmat1@mat
	m2 = scmat2@mat

	miss_nm1 = setdiff(scmat2@genes, scmat1@genes)
	miss_nm2 = setdiff(scmat1@genes, scmat2@genes)

	if (length(miss_nm1) > 0 | length(miss_nm2) > 0) {
	  message(sprintf("will add missing genes (%d to scmat1, %d to scmat2)", length(miss_nm1), length(miss_nm2)))
	}
	m1 = rBind(m1,
					sparseMatrix(i=c(), j=c(),
									dims=c(length(miss_nm1), ncol(m1)),
									dimnames=list(miss_nm1, colnames(m1) )))
	m2 = rBind(m2,
					sparseMatrix(i=c(), j=c(),
									dims=c(length(miss_nm2), ncol(m2)),
									dimnames=list(miss_nm2, colnames(m2) )))
	scmat1@mat = cbind(m1,m2[rownames(m1),])
	scmat1@genes = rownames(m1)
	scmat1@cells= colnames(scmat1@mat)
	scmat1@ncells = ncol(scmat1@mat)
	scmat1@ngenes = nrow(scmat1@mat)
	rownames(scmat1@cell_metadata) = scmat1@cells

	return(scmat1)
}

#'
#' Extract sub-matrix. This return a matrix object on a subset of the genes and cells.
#'
#' @param scmat A tgScMat object.
#' @param genes Genes range. NULL implies the entire gene set
#' @param cells Cells range. NULL implies the entire cell set
#'
#' @export
scm_sub_mat = function(scmat, genes=NULL, cells=NULL)
{
	if (class(scmat)[1] != "tgScMat") {
		stop("invalid scmat in scm_sub_mat")
	}
	if (is.null(genes)) {
		genes = scmat@genes
	}
	if (is.null(cells)) {
		cells = scmat@cells
	}

	if(length(cells)<2) {
		stop("At least 2 cells must be selected")
	}
	if(length(genes)<2) {
		stop("At least 2 genes must be selected")
	}
	meta = scmat@cell_metadata[cells,]
	# remove unused levels
	meta[] <- lapply(meta, function(x) if(is.factor(x)) factor(x) else x)
	return(tgScMat(mat = scmat@mat[genes, cells,drop=F],
						stat_type= scmat@stat_type,
						cell_metadata= meta))
	#batch_map = scmat@batch_map[cells],

}

#' Generate a new matrix object after removing cells without enough umis
#'
#' The currently ignored cells are still going to be ignored, and small cells
#' are goign to be added to them
#'
#' @param new_mat_id id of matrix in scdb
#' @param mat_id existing matrix
#' @param min_umis minimum number of umi per cell
#'
#' @export

mcell_mat_ignore_small_cells = function(new_mat_id, mat_id, min_umis)
{
	if(is.null(scdb_mat(mat_id))) {
		stop("mcell_mat_ignore_small_cells called with missing mat id, ", mat_id)
	}
	mat = scdb_mat(mat_id)
	csize = colSums(mat@mat)
	small_c = names(which(csize < min_umis))
	mcell_mat_ignore_cells(new_mat_id, mat_id, union(mat@ignore_cells, small_c))
}

#' Generate a new matrix object with a given ignore cell list
#'
#' @param new_mat_id id of matrix in scdb
#' @param mat_id existing matrix
#' @param ig_cells cells names to ignore
#' @param reverse set this to true to focus on cells instead of ignore them. False by default
#'
#' @export

mcell_mat_ignore_cells = function(new_mat_id, mat_id, ig_cells, reverse=F)
{
	if(is.null(scdb_mat(mat_id))) {
		stop("mcell_mat_ignore_cells called with missing mat id, ", mat_id)
	}
	mat = scdb_mat(mat_id)

	new_mat = scm_ignore_cells(mat, ig_cells, reverse)
	scdb_add_mat(new_mat_id, new_mat)
}


#' Generate a new matrix object with a given ignore gene list
#'
#' @param new_mat_id id of matrix in scdb
#' @param mat_id existing matrix
#' @param ig_gene gene names to ignore
#' @param reverse set this to true to focus on genes instead of ignore them. False by default
#'
#' @export

mcell_mat_ignore_genes = function(new_mat_id, mat_id, ig_genes, reverse=F)
{
	if(is.null(scdb_mat(mat_id))) {
		stop("mcell_mat_ignore_genes called with missing mat id, ", mat_id)
	}
	mat = scdb_mat(mat_id)

	new_mat = scm_ignore_genes(mat, ig_genes, reverse)
	scdb_add_mat(new_mat_id, new_mat)
}

#' Set ignored (i.e. blacklisted) cells
#'
#' Given a list of cells to ignore, this will cancel any previous policy for blacklisting and remove the given cells to the ignore_mat. Downstream algorithm usually ignore these cells altogether, for any purpose including normalization. However, ignored cells can be accessed and analyzed seperately for validation/tests or when they represent some relevant biology (e.g. cell cycle)
#'
#' @param scmat the matrix object
#' @param ig_cell a list of cell names to ignore
#' @param reverse false by default, if this is true the set of cells to ingore is the complement of the given list
#'
#' @export

scm_ignore_cells = function(scmat, ig_cells, reverse=FALSE)
{
	if(is.null(ig_cells) | length(ig_cells) == 0) {
		ig_cells = vector(l=0)
	}
	if(length(scmat@ignore_cells) > 0) {
		scmat@mat = cbind(scmat@mat, scmat@ignore_cmat)
		if(length(scmat@ignore_genes) > 0) {
			scmat@ignore_gmat = cbind(scmat@ignore_gmat, scmat@ignore_gcmat)
		}
	}
	if(!reverse) {
		good_cells = setdiff(colnames(scmat@mat), ig_cells)
		scmat@ignore_cells = intersect(colnames(scmat@mat), ig_cells)
	} else {
		scmat@ignore_cells = setdiff(colnames(scmat@mat), ig_cells)
		good_cells = intersect(colnames(scmat@mat), ig_cells)
		if(length(good_cells) != length(ig_cells)) {
			stop("Some cells to focus on (ignore_cell, reverse=T), are missing from the current matrix, check your list. len(intersect) = ", length(good_cells), " len(ig_cells) = ", length(ig_cells))
		}
	}
	scmat@cells = good_cells
	scmat@ncells = length(good_cells)
	scmat@ignore_cmat = scmat@mat[,scmat@ignore_cells]
	if(length(scmat@ignore_genes) > 0) {
		scmat@ignore_gcmat = scmat@ignore_gmat[,scmat@ignore_cells]
		scmat@ignore_gmat = scmat@ignore_gmat[,good_cells]
	}
	scmat@mat = scmat@mat[,good_cells]
	return(scmat)
}

#' Set ignored (i.e. blacklisted) genes
#'
#' Given a list of genes to ignore, this will cancel any previous policy for blacklisting and remove the given genes to the ignore_mat. Downstream algorithm usually ignore these genes altogether, for any purpose including normalization. However, ignored genes can be accessed and analyzed seperately for validation/tests or when they represent some relevant biology (e.g. cell cycle)
#'
#' @param scmat the matrix object
#' @param ig_gene a list of gene names to ignore
#' @param reverse false by default, if this is true the set of genes to ingore is the complement of the given list
#'
#' @export

scm_ignore_genes = function(scmat, ig_genes, reverse=FALSE)
{
	if(length(scmat@ignore_genes) > 0) {
		scmat@mat = rbind(scmat@mat, scmat@ignore_gmat)
		if(length(scmat@ignore_cells) > 0) {
			scmat@ignore_cmat = rbind(scmat@ignore_cmat, scmat@ignore_gcmat)
		}
	}
	if(!reverse) {
		good_genes = setdiff(rownames(scmat@mat), ig_genes)
		scmat@ignore_genes = intersect(rownames(scmat@mat), ig_genes)
	} else {
		scmat@ignore_genes = setdiff(rownames(scmat@mat), ig_genes)
		good_genes = intersect(rownames(scmat@mat), ig_genes)
		if(length(good_genes) != length(ig_genes)) {
			stop("Some genes to focus on (ignore_cell, reverse=T), are missing from the current matrix, check your list. len(intersect) = ", length(good_genes), " len(ig_genes) = ", length(ig_genes))
		}
	}
	scmat@genes = good_genes
	scmat@ignore_gmat = scmat@mat[scmat@ignore_genes,]
	if(length(scmat@ignore_cells) > 0) {
		scmat@ignore_gcmat = scmat@ignore_cmat[scmat@ignore_genes,]
		scmat@ignore_cmat = scmat@ignore_cmat[good_genes,]
	}
	scmat@mat = scmat@mat[good_genes,]
	return(scmat)
}

#' Determine recommended downsampling depth.
#'
#' If the parameter scm_n_downsamp is set, it will be returned. If it is null, we return the 5th percentile of the umi count or 750, depending which is larger. However, if the median umi count is lower than this, we use it as the depth. In usch "bad" cases, it may be a good idea to tune the parameter after looking at the umi distribution and resolving possibilities for empty wells/barcodes
#'
#' @param mat_id
#'
#' @export
scm_which_downsamp_n = function(scmat)
{
	downsample_n = get_param("scm_n_downsamp_gstat")
	if(is.null(downsample_n)) {
   	downsample_n = min(round(quantile(colSums(scmat@mat), 0.5)),
                  max(750, round(quantile(colSums(scmat@mat), 0.05))))
  	}
	return(downsample_n)
}

setMethod(
	"show",
	signature = "tgScMat",
	definition =
		function(object) {
			cat("An object of class ", class(object), ",", sep = "")
			cat(" stat type ", object@stat_type, ".\n", sep = "")
			cat(length(object@cells), " cells by ", nrow(object@mat), " genes. median cell content ", median(Matrix::colSums(object@mat)), ".\n", sep = "")
			invisible(NULL)
		}
)

#' Export a mat object (umi matrix and metadata table) to a SingleCellExperiment object
#'
#' @param mat_id id of the mat object to export
#'
#' @return SingleCellExpriment instance
#' @export
#'
scm_export_mat_to_sce = function(mat_id, add_log_counts=F, scale_to=1e4)
{
    mat = scdb_mat(mat_id)
    if (is.null(mat)) {
        error(sprintf("In scm_export_mat_to_seurat, could not find metacell mat %s", mat_id))
    }

    # Adjusted from Seurat Convert functions
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = mat@mat, logcounts = log(1+ scale_to * t(t(mat@mat) / colSums(mat@mat)))))

    SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(mat@cell_metadata[mat@cells, ])

    sce
}

#' Import a umi count matrix with metadata per cell from a SingleCellExperiment objectto a scmat object  to a SingleCellExperiment object
#'
#' @param sce Input SingleCellExperiment to export
#' @param counts_slot name of the umi matrix slot
#'
#' @return tgScMat
#' @export
#'
scm_import_sce_to_mat = function(sce, counts_slot = "counts")
{
    # Adjusted from Seurat Convert functions
    umis <- tryCatch(
        expr = SummarizedExperiment::assay(sce, counts_slot),
        error = function(e) {
            stop(paste0("No data in provided assay - ", counts_slot))
        }
    )

    md <- as.data.frame(SummarizedExperiment::colData(sce))

    scm_new_matrix(umis, md)
}
