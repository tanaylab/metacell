#' Gene sets interface
#'
#' Gene sets are simply container for multiple sets of genes IDs. This is in use when extracting features for meta-cell analysis, when deriving modules of correlated genes, or when importing annotated genes or transferring them between analyses (e.g. cell cycle).
#'
#' @slot description textual description of the gene set collection
#' @slot set_names possible annotations of gene sets
#' @slot gene_set set per gene
#'
#' @import Matrix
#'
#' @export tgGeneSets
#' @exportClass tgGeneSets
tgGeneSets <- setClass(
   "tgGeneSets",
	slots = c(
	  description = "character",
	  set_names = "vector",
	  gene_set = "vector")
)

#' tgGeneSets public constructor
#'
#' This constructs a tgGene sets object, nothing fancy
#'
#' @param sets a named vector assigning set id to genes
#' @param desc an optional string describing the gene set
#' @export

setMethod(
  "initialize",
  signature = "tgGeneSets",
  definition =
    function(.Object, sets, desc=NULL) {
      .Object@description= ifelse(is.null(desc), "NS", desc)
		set_names = unique(sets)
		.Object@set_names = set_names
		.Object@gene_set = sets
      return(.Object)
    }
)

#' Generating a new gset
#
#' @param sets a vector where names are genes and values are set ids
#' @param desc tetual description of the gset
#' @export
gset_new_gset = function(sets, desc)
{
	return(tgGeneSets(sets, desc))
}

#' Add specific genes to the set
#
#' @param gset gene set object
#' @param genes genes to add
#' @param subset_id subset id to which the new genes will be assigned
#' @export
gset_add_genes = function(gset, genes, subset_id)
{
	gset@gene_set[genes] = subset_id
	return(gset)
}

#' Import gene set for a text table
#
#' @param fn file name to read from
#' @param desc tetual description of the gset, if null the desc will be the file name
#' @export
gset_import_table = function(fn, desc=NULL)
{
	sets = read.table(fn, h=T, sep="\t")
	if(is.null(desc)) {
		desc = fn
	}
	if("gene" %in% colnames(sets) & "set" %in% colnames(sets)) {
		gs = sets$set
		names(gs) = set$gene
		return(tgGeneSets(gs, desc))
	} else {
		stop("cannot initialize gset from ", fn, " file must be a tab delim table with a header including the fields gene and set\n")
	}
}

#' Exprt gene set to a table
#
#' @param gset gene set
#' @param fn file name to save to
#' @export
gset_write_table = function(gset, fn)
{
	write.table(data.frame(gene = names(gset@gene_set), set = gset@gene_set), fn, quote=F, sep="\t")
}

#' Get genes of one set
#
#' @param gset a gene set object
#' @param id the set id
#' @export
gset_get_genes = function(gset, id)
{
	return(names(gset@gene_set)[gset@gene_set == id])
}

#' Generate a new gene set from an existing one, filtered by a list of genes
#
#' @param gset a gene set object
#' @param genes names of genes to filter by
#' @param inverse if true, omit the genes in the genes parameter
#' @param desc description of the new gene set
#' @export
gset_new_restrict_nms = function(gset, genes, inverse=F, desc)
{
	sets = gset@gene_set
	if(inverse) {
		rsets = sets[setdiff(names(sets),genes)]
	} else {
		rsets = sets[intersect(names(sets),genes)]
	}

	return(tgGeneSets(rsets, desc))
}

#' Generate a new gene set from an existing one, filtered by a list of genes
#
#' @param gset a gene set object
#' @param filt_gset gset of genes to filter by
#' @param inverse if true, omit the genes in the genes parameter
#' @param desc description of the new gene set
#' @export
gset_new_restrict_gset = function(gset, filt_gset, inverse=F, desc)
{
	sets = gset@gene_set
	genes = names(filt_gset@gene_set)
	if(inverse) {
		rsets = sets[setdiff(names(sets),genes)]
	} else {
		rsets = sets[intersect(names(sets),genes)]
	}

	return(tgGeneSets(rsets, desc))
}

#' extract umi matrix for the genes in the set, possibly donwsampling
#'
#' @param gset_id gene set in scdb
#' @param mat_id mat in scdb
#' @param downsamp if this is true the returned matrix is downsampled
#'
#' @export
gset_get_feat_mat = function(gset_id, mat_id, downsamp = F)
{
	gset = scdb_gset(gset_id)
	if(is.null(gset)) {
		stop("MC-ERR non existing gset in gset_get_feat_mat, id ", gset_id)
	}
	mat = scdb_mat(mat_id)
	if(is.null(mat)) {
		stop("MC-ERR non existing mat in gset_get_geat_mat, id ", mat_id)
	}

	if(downsamp) {
		downsample_n = scm_which_downsamp_n(mat)
		message("will downsample the matrix, N= ", downsample_n, " (and yes - this should have been chached")
		umis = scm_downsamp(mat@mat, downsample_n)
	} else {
		umis = mat@mat
	}
	nms = names(gset@gene_set)
	nms = intersect(nms, rownames(umis))
	if(length(nms) == 0) {
		stop("get feature matrix with zero overlap of gene names with gset_id ", gset_id, " mat ", mat_id)
	}
	feat_ds = umis[nms,]
	return(feat_ds)
}
