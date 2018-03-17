#' Meta cell cover
#'
#' Representing a meta cell cove of a given cell graph (or more generally of a scRNA data matrix)
#'
#' @slot mc assignment of cells to metacell id
#' @slot cell_names names of cells (all other objects use running integres relating to these names)
#' @slot mc_fp a matrix showing for each gene (row) the relative enrichment of umis
#' @slot e_gc a matrix showing for each gene (row) the mean umi count (without normalizing for cell size)
#' @slot cov_gc a matrix showing for each gene (row) the fraction of cells in the meta cell that are non zero for the umi
#' @slot n_bc a matrix detemrining for each batch (row) the meta cell breakdown
#' @slot annots names of metacells (oridnal numbers by default)
#' @slot colors colors of metacells
#' @slot color_key data.frame defining color per markers genes
#'
#' @import Matrix
#'
#' @export tgMCCov
#' @exportClass tgMCCov
tgMCCov <- setClass(
   "tgMCCov",
	slots = c(
	  mc = "vector",
	  outliers = "vector",
	  cell_names = "vector",
	  mc_fp = "matrix",
	  e_gc = "matrix",
	  cov_gc = "matrix",
	  n_bc = "matrix",
	  annots = "vector",
	  colors = "vector",
	  color_key = "data.frame")
)

#' Construct a meta cell object 
#'
#' This constructs a meta cell cover object. It gets an MC assignment (cell->MC_ID), and a matrix, and call standard api of this class to compute the footprints.
#'
#' @param mc assignment of metacell id to cell
#' @param scmat a single cell RNA matrix object
#' @export

setMethod(
  "initialize",
  signature = "tgMCCov",
  definition =
    function(.Object, mc, outliers = c(), scmat) {
		all_cells = colnames(scmat@mat)
		.Object@mc = mc
		.Object@outliers = outliers
		.Object@cell_names = all_cells
		mc_cells = names(mc)
		if(length(intersect(outliers, mc_cells)) > 0) {
			stop("MC-ERR non zero intersect of outliers and mc assigned cells")
		}
		if(length(all_cells) != (length(mc_cells)+length(outliers))
		| length(setdiff(all_cells, c(mc_cells, outliers))) > 0) {
			stop("matrix cell names and mc assignments+outliers differ in tgMCCov initialization")
		}
		max_clust = max(mc)
		.Object@colors = rep("white", n=max_clust)
		.Object@annots= 1:max_clust
		.Object = mc_update_stats(.Object, scmat)
		.Object@color_key = data.frame(group=c(), gene=c(), color=c(), priority=c(), T_fold=c())
      return(.Object)
    }
)

#' Generate a new metacell in scdb
#'
#' This constructs a meta cell cover object and puts it into scdb. It gets an MC assignment (cell->MC_ID), and a matrix, and call standard api of this class to compute the footprints.
#'
#' @param mc_id id of scdb meta cell object ot be added
#' @param mc assignment of metacell id to cell
#' @param scmat a single cell RNA matrix object
#' @export
mcell_new_mc = function(mc_id, mc, outliers, scmat)
{
	scdb_add_mc(mc_id, tgMCCov(mc, outliers, scmat))
}

#' Compute stats over metacell and update the object
#'
#' This compute statistics per metacell given an sc matirx.
#'
#' @param mc a metacell object
#' @param scmat a matrix object to draw statisics from
#' @export
mc_update_stats = function(mc, scmat)
{
	us = scmat@mat[, names(mc@mc)]
	message("add batch counts")
	mc@n_bc = mc_compute_n_bc(mc, scmat)
	message("compute footprints")
	mc@mc_fp = mc_compute_fp(mc, us)
	message("compute absolute ps")
	mc@e_gc= mc_compute_e_gc(mc, us)
	message("compute coverage ps")
	mc@cov_gc = mc_compute_cov_gc(mc, us)
	return(mc)
}

#' Move all cels from specific metacells to the outliers
#'
#' This can be used to treat doublets metacells or other artifacts.
#'
#' @param mc mc object
#' @param mc_id metacell to eliminate
#'
#' @export
mc_set_outlier_mc = function(mc, mc_id)
{
	N = ncol(mc@mc_fp)
	id_map = rep(-1,N)
	id_map[(1:N)[-mc_id]] = 1:(N-1)
	ocells = names(mc@mc)[which(mc@mc==mc_id)]
	mc@outliers = c(mc@outliers, ocells)
	gcells = setdiff(mc@cell_names, mc@outliers)
	mc@mc = mc@mc[gcells]
	mc@mc = id_map[mc@mc]
	names(mc@mc) = gcells
	drop_f = colnames(mc@mc_fp) != mc_id
	mc@mc_fp = mc@mc_fp[,drop_f]
	colnames(mc@mc_fp) = 1:(N-1)
	mc@e_gc = mc@e_gc[,drop_f]
	mc@cov_gc = mc@cov_gc[,drop_f]
	colnames(mc@cov_gc) = 1:(N-1)
	mc@n_bc = mc@n_bc[,drop_f]
	colnames(mc@n_bc) = 1:(N-1)
	mc@annots = mc@annots[drop_f]
	names(mc@annots) = 1:(N-1)
	mc@colors= mc@colors[drop_f]
	names(mc@colors) = 1:(N-1)
	return(mc)
}


#' Compute metacell gene footprint 
#'
#' The footprint is defined as the size-normalized geometric mena of the number of umis per metacells, dvided by the median over all metacells.
#'
#' @param mc a metacell object
#' @param us umi matrix
#' @export
mc_compute_fp = function(mc, us)
{
	f_g_cov = rowSums(us) > 10

	clust_geomean = .row_stats_by_factor(us[f_g_cov,], 
									mc@mc, 
									function(y) {exp(rowMeans(log(1+y)))-1})

	clust_meansize = tapply(colSums(us), mc@mc, mean)
	ideal_cell_size = pmin(1000, median(clust_meansize))
	g_fp = t(ideal_cell_size*t(clust_geomean)/as.vector(clust_meansize))
	#normalize each gene
	fp_reg = 0.1   
	#0.1 is defined here because 0.1*mean_num_of_cells_in_cluster 
	#is epxected to be 3-7, which means that we regulairze 
	#umicount in the cluster by 3-7.
	g_fp_n = (fp_reg+g_fp)/apply(fp_reg+g_fp, 1, median)

	return(g_fp_n)
}

#' Compute metacell absolute mean umi per cell
#'
#' This compute the genometric mean of the number of umis per cells for each metacell
#'
#' @param mc a metacell object
#' @param us umi matrix
#' @export
mc_compute_e_gc= function(mc, us)
{
	f_g_cov = rowSums(us) > 10

	e_gc = .row_stats_by_factor(us[f_g_cov,], 
									mc@mc, 
									function(y) {exp(rowMeans(log(1+y)))-1})

	return(e_gc)
}

#' Compute fraction of non zero expressing cells per gene and mc
#'
#'
#' @param mc a metacell object
#' @param us umi matrix
#' @export
mc_compute_cov_gc= function(mc, us)
{
	cov_gc = .row_stats_by_factor(us > 0, 
										fact = mc@mc, 
										rowFunction = rowMeans)
	return(cov_gc)
}

#' Compute distribution of cells over batches and metacell
#'
#' @param mc a metacell object
#' @param scmat scamt object (for the metadata)
#' @export
mc_compute_n_bc= function(mc, scmat)
{
	tb = table(scmat@cell_metadata[names(mc@mc), "amp_batch_id"], mc@mc)
	n_bc = matrix(tb, ncol=dim(tb)[2])
	rownames(n_bc) = dimnames(tb)[[1]]
	return(n_bc)
}

#' Update metacell annotation
#'
#' @param mc a metacell object
#' @param annot annotation vectors
#' @export
mcell_mc_add_annot = function(mc_id, annots)
{
	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
		stop("MC-ERR: missing mc_id when updating annotation id = ", mc_id)
	}
	mc@annots = annots
	scdb_add_mc(mc_id, mc)
}

#' Update metacells  colors 
#'
#' @param mc a metacell object
#' @param colors a vector of colors per metacell
#' @export
mcell_mc_add_color= function(mc_id, colors)
{
	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
		stop("MC-ERR: missing mc_id when updating annotation id = ", mc_id)
	}
	mc@colors = colors 
	scdb_add_mc(mc_id, mc)
}

#' Reorder metacell data given defined order
#'
#' @param mc metacell object
#' @param ord new order metacells
#' 
#' @export

mc_reorder = function(mc, ord)
{
	if(length(ord) != ncol(mc@mc_fp)) {
		stop("MC-ERR: reordering metacells with an order vector shorter than the number of MC in the object")
	}
	ord_map = rep(NA, ncol(mc@mc_fp))
	ord_map[ord] = 1:ncol(mc@mc_fp)
	if(sum(is.na(ord_map)) > 0) {
		stop("MC-ERR: reordering metacells with an order vector lacking some of the ids - generate a new subset metacell object if this is what you wanted to do")
	}
	nms = names(mc@mc)
	mc@mc = ord_map[mc@mc]
	names(mc@mc) = nms
	mc@mc_fp = mc@mc_fp[,ord]
	colnames(mc@mc_fp) = 1:length(ord)
	mc@e_gc = mc@e_gc[,ord]
	colnames(mc@e_gc) = 1:length(ord)
	mc@cov_gc = mc@cov_gc[,ord]
	colnames(mc@cov_gc) = 1:length(ord)
	bnms = rownames(mc@n_bc)
	mc@n_bc = matrix(mc@n_bc[,ord], ncol=length(ord))	#as.matrix to avoid casting on 1 batch data
	rownames(mc@n_bc) = bnms
	colnames(mc@n_bc) = 1:length(ord)
	if(!is.null(mc@annots)) {
		mc@annots = mc@annots[ord]
		names(mc@annots) = 1:length(mc@annots)
	}
	if(!is.null(mc@colors)) {
		mc@colors = mc@colors[ord]
		names(mc@colors) = 1:length(mc@colors)
	}
	return(mc)	
}

#' Reorder metacells using hierarchical clustering 
#'
#' MEtacells are reorder according to their footprint similarity based on hclust and the reordering using two select antagonistic markers (that can be selected automatically)
#'
#' @param mc_id id of metacell object
#' @param gene_left gene for reordering toward the left side (null by default)
#' @param gene_right gene for reordering toward the right side (null by default)
#' @param gset_blist_id id of gene set to blacklist while ordering (e.g. cell cycle genes)
#'
#'

mcell_mc_reorder_hc = function(mc_id, 
							gene_left = NULL, gene_right = NULL, 
							gset_blist_id = NULL)
{
	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
		stop("MC-ERR: missing mc_id when updating annotation id = ", mc_id)
	}
	gene_folds = mc@mc_fp

	marks = rownames(gene_folds)

	if(!is.null(gset_blist_id)) {
		blist = scdb_gset(gset_blist_id)
		if(is.null(blist)) {
			stop("MC-ERR unknown geneset id ", gset_blist_id, " when blacklisting markers for ordering metacells")
			marks = setdiff(marks, names(blist@gene_set))
		}
	}
	max_gcov = apply(mc@cov_gc[marks,],1,max)
	cov_marks = rownames(mc@cov_gc)[max_gcov > 0.25]
	if(length(intersect(marks, cov_marks)) > ncol(gene_folds)) {
		marks = intersect(marks, cov_marks)
	}

	gene_folds = gene_folds[marks,]

	good_marks = unique(as.vector(unlist(
			apply(gene_folds,
				2,
				function(x)  {
				   names(head(sort(-x[x>0.5]),n=10)) })
		     )))
	
	if(is.null(good_marks) | length(good_marks) < 4) {
		good_marks= rownames(gene_folds)
	}
	feat = log2(gene_folds[good_marks,])

	hc = hclust(dist(cor(feat)), "ward.D2")

	if(!is.null(gene_right) & !is.null(gene_left)) {
		if(!gene_right %in% marks | !gene_left %in% marks) {
			stop("genes for order polarization ", gene_right, " ", gene_left, " are not in markers")
		}
	}
	if(is.null(gene_right) | is.null(gene_left)) {
		g_ncover = apply(feat > 1, 1, sum)
		main_mark = names(g_ncover)[which.max(g_ncover)]
		f = feat[main_mark,] < 0.25
		if(sum(f) > 0.2*ncol(feat)) {
			g_score = apply(feat[,f]>1, 1, sum)
		} else {
			g_score = -apply(feat, 1, cor, feat[main_mark,])
		}
		second_mark = names(g_score)[which.max(g_score)]
		message("reorder on ", main_mark, " vs ", second_mark)
	}

	d = reorder(as.dendrogram(hc), 
				feat[gene_right,]-feat[gene_left,], 
				agglo.FUN=mean)
	hc2 = as.hclust(d)

	scdb_add_mc(mc_id, mc_reorder(mc, hc2$order))
}

#' TEst if a graph object cover all cells in the mc
#'
#' mc_id mc object
#' graph_id graph object
#'
#' @export
#'
mcell_mc_match_graph = function(mc_id, graph_id) 
{
	mc = scdb_mc(mc_id)
	graph = scdb_cgraph(graph_id)
	cov = intersect(graph@cell_names, mc@cell_names)
	return(length(cov) == length(mc@cell_names)) 
}