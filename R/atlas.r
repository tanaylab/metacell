#' wrap up an atlas object from mc, mc2d and matrix ids
#'
#' This is not doing much more than generating a list with the relevant object names bundeled. To be enhnaced at some stage.
#'
#' @param mat_id id of metacell object ina scdb
#' @param mc_id id of metacell object ina scdb
#' @param gset_id features defining the atlas (to be sued for determing projection)
#' @param mc2d_id projection object id, to define the atlas 2D layout
#'
#' @export

#' @export tgMCAtlas

tgMCAtlas <- setClass(
   "tgMCAtlas",
	slots = c(
	  atlas_name = "character",
	  mat_id = "character",
	  gene_naming_type = "character",
	  mc_id = "character",
	  gset_id = "character",
	  cgraph_id = "character",
	  edges_Knn = "integer",
	  cell_cell_mean_cor = "vector",
	  cell_mc_cor = "vector",
	  cell_adj_type_n = "matrix",
	  color_type_id = "vector",
	  type_color = "vector",
	  type_name = "vector")
)

#' Construct a meta cell reference atlas
#'
#'
#' @param mat_id umi matrix of the reference
#' @param naming_type gene naming type 
#' @param mc_id metacell object id 
#' @param graph_id KNN adjacencies for the atlas
#' @param gset_id gene set id for defining features
#' @param type_annot data frame with type_id, type_name, type_color fields
#' @export

setMethod(
  "initialize",
  signature = "tgMCAtlas",
  definition =
    function(.Object, mat_id, naming_type,
							 mc_id, graph_id, gset_id, 
							 type_annot) {
		.Object@atlas_name = mat_id
		.Object@mat_id = mat_id
		.Object@mc_id = mc_id
		.Object@gset_id = gset_id 
		.Object@cgraph_id = graph_id 
		.Object@gene_naming_type = naming_type
		color_type_id = type_annot$type_id
		names(color_type_id) = toupper(type_annot$type_color)
		type_name = as.character(type_annot$type_name)
		names(type_name) = as.character(type_annot$type_id)
		type_color = as.character(type_annot$type_color)
		names(type_color) = as.character(type_annot$type_id)
		.Object@color_type_id = color_type_id
		.Object@type_name = type_name 
		.Object@type_color = type_color 
		mc = scdb_mc(mc_id)
		if(is.null(mc)) {
			stop("MC-ERR unkown mc_id ", mc_id, " when building atlas")
		}
		mat = scdb_mat(mat_id)
		if(is.null(mat)) {
			stop("MC-ERR unkown mat_id ", mat_id, " when building atlas")
		}
		graph = scdb_cgraph(graph_id)
		if(is.null(graph)) {
			stop("MC-ERR unkown graph_id ", graph_id, " when building atlas")
		}
		.Object@edges_Knn = max(table(graph@edges$mc1))

		adj_type_n = mcatlas_comp_adj_type_n(mc, graph, color_type_id)
		.Object@cell_adj_type_n = as.matrix(adj_type_n)
		.Object = mcatlas_update_internal_cors(.Object)

      return(.Object)
    }
)

#' Generate a new atlas in scdb
#'
#' This constructs a meta cell atlas object 
#'
#' @param mat_id id matrix object
#' @param naming_type gene name type
#' @param mc_id id of metacell object
#' @param graph_id id of cgraph object
#' @param gset_id id of feature gene set
#' @param edges_Knn the knn parameter for neighbor type distribution
#' @param type_annot data frame with type_id, type_name, type_color fields
#' @export
mcell_new_mcatlas = function(atlas_id, mat_id, naming_type, mc_id, graph_id, gset_id, type_annot)
{
	scdb_add_mcatlas(atlas_id, 
					tgMCAtlas(mat_id, 
								naming_type=naming_type, 
								mc_id=mc_id, graph_id=graph_id, gset_id = gset_id, 
								type_annot = type_annot))
}

#' Compute the summary matrix of cell type adjacencies 
#'
#'
#' @param mc metacell object for defining colors/type per cell
#' @param cgraph cell-cell knn graph
#' @param color_ype_id named vector mapping colors to IDs
#'
#' @return summary matrix with cells in rows, types in columns
#' @export

mcatlas_comp_adj_type_n = function(mc, graph, color_type_id)
{
	gr = graph@edges
	col_tab= toupper(mc@colors[mc@mc[levels(gr$mc1)]])
	gr$col1 = col_tab[gr$mc1]
	gr$col2 = col_tab[gr$mc2]

	n1 = gr %>% group_by(mc1, col2) %>% summarise(tot=n()) %>% spread(col2, value="tot", fill=0)
	n2 = gr %>% group_by(mc2, col1) %>% summarise(tot=n()) %>% spread(col1, value="tot", fill=0)

	n1 = as.data.frame(n1)
	n2 = as.data.frame(n2)

	rownames(n1) = n1[,1]
	rownames(n2) = n2[,1]
	cells = names(mc@mc)
	act_types = intersect(colnames(n1), names(color_type_id))
	act_types = intersect(colnames(n2), act_types)
	n1 = n1[cells, act_types]
	n2 = n2[cells, act_types]
	n1[is.na(n1)] = 0
	n2[is.na(n2)] = 0
	adj_n = n1+n2
	colnames(adj_n) = color_type_id[colnames(adj_n)]
	return(adj_n)
}

mcatlas_update_internal_cors = function(atlas)
{
	mc = scdb_mc(atlas@mc_id)
	mat = scdb_mat(atlas@mat_id)
	cgraph = scdb_cgraph(atlas@cgraph_id)
	gset = scdb_gset(atlas@gset_id)

	gnms = names(gset@gene_set)
	gnms = intersect(gnms, rownames(mc@e_gc))
	mc_feats = log2(mc@e_gc[gnms,]+1e-5)

	dsamp_n = quantile(colSums(mat@mat),0.1)
	ds_umis = gset_get_feat_mat(gset_id=atlas@gset_id, mat_id=atlas@mat_id, 
							downsamp = T, add_non_dsamp=T,
							downsample_n = dsamp_n)
	k_nonz_exp = get_param("scm_k_nonz_exp", "metacell")
	sc_feats = log2(1+k_nonz_exp*ds_umis)

	message("Will compute Knn correlation thresholds for all atlas cells")
	sc_cor = tgs_cor_knn(x=as.matrix(sc_feats), y=as.matrix(sc_feats), knn=atlas@edges_Knn)
	message("Will compute best mc correlation to atlas MC")
	mc_cor = tgs_cor_knn(as.matrix(sc_feats), as.matrix(mc_feats), knn=1)

	atlas@cell_mc_cor = mc_cor$cor
	names(atlas@cell_mc_cor) = as.character(mc_cor$col1)

	mean_c = tapply(sc_cor$cor, sc_cor$col1, mean)	
	atlas@cell_cell_mean_cor = mean_c

	return(atlas)	
}

#' Compute association of query umi matrix and atlas
#'
#'
#' @param atlas and atlas object
#' @param qmat query matrix
#' @param qmat_naming_type gene naming type of query
#'
#' @return a query processed list 
#' @export

mcatlas_project_query = function(atlas, qmat_id, qmat_naming_type, all_vs_all)
{
	a_mc = scdb_mc(atlas@mc_id)
	#merge matrices
	K = atlas@edges_Knn
	color_type_id = atlas@color_type_id
	names(color_type_id) = toupper(names(color_type_id))
	qmat = scdb_mat(qmat_id)
	if(is.null(qmat)) {
		stop("Q matrix id ", qmat_id, " not found in scdb")
	}
	qmat_dsamp_n = quantile(colSums(qmat@mat),0.1)

	amat = scdb_mat(atlas@mat_id)
	amat_dsamp_n = quantile(colSums(amat@mat),0.1)

	dsamp_n = floor(min(amat_dsamp_n, qmat_dsamp_n))

	gene_name_src_targ = NULL
	if(atlas@gene_naming_type != qmat_naming_type) {
		gene_name_src_targ = c(qmat_naming_type, atlas@gene_naming_type)
	}

	q_feats = gset_get_feat_mat(gset_id=atlas@gset_id, mat_id=qmat_id, 
							downsamp =T, add_non_dsamp=T,
							downsample_n = dsamp_n, 
							gene_names_src_targ=gene_name_src_targ)

	a_feats = gset_get_feat_mat(gset_id=atlas@gset_id, mat_id=atlas@mat_id, 
							downsamp = T, add_non_dsamp=T,
							downsample_n = dsamp_n)
	
	if(nrow(q_feats) < nrow(a_feats)*0.3) {
		stop("less than 30% of the atlas features are mapped in query genes - improve naming type conversion table?")
	}
	a_feats = a_feats[rownames(q_feats),]

	feat = cbind(q_feats, a_feats)
	
	k_nonz_exp = get_param("scm_k_nonz_exp", "metacell")

	if(all_vs_all) {
		feat = log2(1+k_nonz_exp*as.matrix(feat))

		message("will build balanced knn graph on ", ncol(feat), " cells and ", nrow(feat), " genes, this can be a bit heavy for >20,000 cells")
		k_alpha = get_param("scm_balance_graph_k_alpha", "metacell")
		k_beta = get_param("scm_balance_graph_k_beta", "metacell")
		k_expand = get_param("scm_balance_graph_k_expand", "metacell")
		gr = tgs_cor_graph(x=feat, knn=K, k_expand=k_expand, k_alpha=k_alpha, k_beta=k_beta)

		query_cells = colnames(q_feats)

		colnames(gr) = c("cell1","cell2","w")

		col_tab= toupper(a_mc@colors[a_mc@mc[levels(gr$cell1)]])
		is_query = levels(gr$cell1) %in% query_cells

		gr$col1 = col_tab[gr$cell1]
		gr$col2 = col_tab[gr$cell2]

		gr$query1 = is_query[gr$cell1]
		gr$query2 = is_query[gr$cell2]

		q_gr1 = gr[gr$query1,]
		q_gr2 = gr[gr$query2,]

		q_nadj = tapply(c(q_gr1$query2,q_gr2$query1), 
												c(q_gr1$cell1,q_gr2$cell2), length)
		q_nadj_q = tapply(c(q_gr1$query2,q_gr2$query1), 
												c(q_gr1$cell1,q_gr2$cell2), sum)

		type_n1 = q_gr1 %>% group_by(cell1, col2) %>% summarise(tot=n()) %>% spread(col2, value="tot", fill=0)
		type_n2 = q_gr2 %>% group_by(cell2, col1) %>% summarise(tot=n()) %>% spread(col1, value="tot", fill=0)

		type_n1 = as.data.frame(type_n1)
		type_n2 = as.data.frame(type_n2)
		rownames(type_n1) = type_n1[,1]
		rownames(type_n2) = type_n2[,1]
		act_types = intersect(colnames(type_n1), names(color_type_id))
		act_types = intersect(colnames(type_n2), act_types)
		type_n1 = type_n1[query_cells, act_types]
		type_n2 = type_n2[query_cells, act_types]
		type_n1[is.na(type_n1)] = 0
		type_n2[is.na(type_n2)] = 0
		q_adj_n = type_n1+type_n2
		colnames(q_adj_n) = color_type_id[colnames(q_adj_n)]

		return(list(q_adj_n = q_adj_n, q_mean_cor = NULL, 
							q_nadj = q_nadj, q_nadj_q = q_nadj_q))
	} else {
		gr = tgs_cor_knn(log2(1+k_nonz_exp*as.matrix(q_feats)),
									log2(1+k_nonz_exp*as.matrix(a_feats)), K)

		query_cells = colnames(q_feats)

		colnames(gr) = c("cell1","cell2","cor", "rank")

		mean_cor = tapply(gr$cor, gr$cell1, mean, na.rm=T)	
		ref_mean_cor = tapply(atlas@cell_cell_mean_cor[gr$cell2], gr$cell1, mean, na.rm=T)

		col_tab= toupper(a_mc@colors[a_mc@mc[levels(gr$cell2)]])
		gr$col2 = col_tab[gr$cell2]

		type_n2 = gr %>% group_by(cell1, col2) %>% summarise(tot=n()) %>% spread(col2, value="tot", fill=0)

		type_n2 = as.data.frame(type_n2)
		rownames(type_n2) = type_n2[,1]
		act_types = intersect(colnames(type_n2), names(color_type_id))
		type_n2 = type_n2[query_cells, act_types]
		type_n2[is.na(type_n2)] = 0
		q_adj_n = type_n2
		colnames(q_adj_n) = color_type_id[toupper(colnames(q_adj_n))]
		return(list(q_adj_n = q_adj_n, q_mean_cor = mean_cor, 
									q_ref_mean_cor = ref_mean_cor,
											q_nadj = NULL, q_nadj_q = NULL))
	}
}

mcatlas_annotate_sc_by_projection = function(atlas, qmat_id, qmat_naming_type, all_vs_all)
{
	res = mcatlas_project_query(atlas, qmat_id, qmat_naming_type, all_vs_all)

	a_mc = scdb_mc(atlas@mc_id)
#find best hit
	no_hits_nms = setdiff(colnames(atlas@cell_adj_type_n), colnames(res$q_adj_n))
	ndiff = length(no_hits_nms)
	if(ndiff > 0) {
		q_adj_n = cbind(res$q_adj_n, matrix(0, nrow=nrow(res$q_adj_n), ncol=ndiff))
		colnames(q_adj_n) = all_nms
	} else {
		q_adj_n = res$q_adj_n
	}
	match = tgs_cor_knn(t(res$q_adj_n), t(atlas@cell_adj_type_n[,names(res$q_adj_n)]),knn=1)

	names(match) = c("q_sc_id", "a_sc_id", "cor", "rank")
	match$annot_col = toupper(a_mc@colors[a_mc@mc[match$a_sc_id]])

#is it good enough hit?	
	qsc_type = atlas@color_type_id[as.character(match$annot_col)]
	names(qsc_type) = match$q_sc_id
	annots = data.frame(cell_id = match$q_sc_id, 
							  color = as.character(match$annot_col),
							  type = qsc_type, 
							  name = atlas@type_name[qsc_type])
	if(!is.null(res$q_mean_cor)) {
		 annots$cor = res$q_mean_cor
		 annots$ref_top_cor = res$q_ref_mean_cor
	} 
	if(!is.null(res$q_nadj_q)) {
		 annots$ref_pure = res$q_nadj_q/res$q_nadj
	}

	return(annots)
}

#' Annotate query metacell with atlas by comparing query cells  to altas cells with a Knn strategy

#' Algorithm short description TBA
#'
#' @param atlas_id id of atlas object in scdb
#' @param qmc_id id of metacell object ina scdb
#' @param qmat_id id of query umi matrix
#' @param qmat_naming_type naming scheme of query matrix/mc
#' @param all_vs_all should a graph combining query and atlas be constructed, or only K-nn of a bipartiate query atlas graph be examined.
#' @param new_qmc_id id of recolored metacell object to save in DB (NULL will supress updating db)
#' @param fig_cmp_dir name of directory to put figures per MC 
#' @param T_cor_gap how much gap between internal atlas mc-mc correlation and query-atlas correlation one is allowing to keep the annotation.
#' @param sc_annot cached sc_annot object can be specified (as returned by mcatlas_annotate_sc_by_projection) - to save time in repeated annotation refinement
#'
#' @export

mcatlas_annotate_mc_by_sc2sc_projection = function(atlas_id, 
					qmc_id, qmat_id, qmat_naming_type, 
					all_vs_all, new_qmc_id, 
					T_cor_gap=1, fig_cmp_dir=NULL, sc_annots=NULL)
{
	atlas = scdb_mcatlas(atlas_id)
	if(is.null(atlas)) {
		stop("cannot find atlas id ", atlas_id, " when callig annotate mc by sc projection")
	}

	fig_cmp_dir = gen_fig_cmp_dir(fig_cmp_dir, qmc_id, atlas@atlas_name)

	if(is.null(sc_annots)) {
		sc_annots = mcatlas_annotate_sc_by_projection(atlas, qmat_id, qmat_naming_type, all_vs_all)
		sc_annots$color = as.character(sc_annots$color)
	} else {
		message("using cached sc annots")
	}

	mc = scdb_mc(qmc_id)

	mc_type_ns = table(mc@mc, sc_annots[names(mc@mc),"type"])
	mc_type_ns_n = mc_type_ns/rowSums(mc_type_ns)
	mc_type = colnames(mc_type_ns_n)[apply(mc_type_ns_n, 1, which.max)]

	cor_gap = sc_annots[names(mc@mc),"cor"]-sc_annots[names(mc@mc), "ref_top_cor"]
	mc_avg_cor = tapply(sc_annots[names(mc@mc),"cor"], mc@mc, mean)
	mc_avg_ref_cor = tapply(sc_annots[names(mc@mc), "ref_top_cor"], mc@mc, mean)
	mc_avg_cor_gap  = mc_avg_cor - mc_avg_ref_cor

	mc@colors = ifelse(mc_avg_cor_gap > -T_cor_gap, as.character(atlas@type_color[mc_type]), "gray")

#plot scatter of correlation hit and internal atlas correlation on self hits
	png(sprintf("%s/proj_cor_vs_atlas_self.png", fig_cmp_dir),w=600,h=600)
	plot(sc_annots[names(mc@mc),"cor"], sc_annots[names(mc@mc), "ref_top_cor"], pch=19, col=sc_annots[names(mc@mc),"color"], cex=0.1)
	points(mc_avg_cor, mc_avg_ref_cor, pch=19, cex=0.8, col=mc@colors)
	abline(a=0, b=1)
	abline(a=T_cor_gap, b=1)
	dev.off()
#plot bars of best proj color on each mc + entropy
	mc_proj_col_p = mc_type_ns_n
	colnames(mc_proj_col_p) = atlas@type_color[colnames(mc_type_ns_n)]
	query_entropy = apply(mc_proj_col_p, 1, 
				function(x) { p =x/sum(x); lp =log2(p+1e-6); return(sum(-p*lp)) })

	png(sprintf("%s/sc2sc_query_color_dist.png", fig_cmp_dir),w=600,h=150+12*nrow(mc_proj_col_p))
	layout(matrix(c(1,2),nrow=1),widths=c(5,2))
	par(mar=c(2,3,2,0))
	atlas_cols = intersect(atlas@type_color,colnames(mc_proj_col_p))
	barplot(t(mc_proj_col_p[,atlas_cols]), col=atlas_cols, horiz=T, las=2)
	par(mar=c(2,0,2,3))
	barplot(query_entropy, col=mc@colors, horiz=T, las=2)
	grid()
	dev.off()
	
	scdb_add_mc(new_qmc_id, mc)
	return(sc_annots)
}


#' Annotate query metacell with atlas by comparing query metacell gene profiles to the atlas MCs gene profiles

#' This will take each MC in the  query MC and find its to correlated
#' metacell in the reference, gemeratomg some figures along the way/
#'
#' @param atlas_id id of atlas object in scdb
#' @param qmc_id id of metacell object ina scdb
#' @param qmat_naming_type naming scheme of query matrix/mc
#' @param new_qmc_id id of recolored metacell object to save in DB (NULL will supress updating db)
#' @param fig_cmp_dir name of directory to put figures per MC 
#' @param q_gset_id query gene set id object (optional) - to restrict features used in comaprison to the intersection of atlas and query gsets
#' @param T_cor_gap how much gap between internal atlas mc-mc correlation and query-atlas correlation one is allowing to keep the annotation.
#' @param nd_color - color for metacell without an annotation
#'
#' @export

mcatlas_annotate_mc_by_mc2mc_projection = function(atlas_id, 
					qmc_id, qmat_naming_type, 
					new_qmc_id, 
					q_gset_id = NULL,
					T_cor_gap=1, nd_color = "lightgray", fig_cmp_dir=NULL)
{
#using only the 
	atlas = scdb_mcatlas(atlas_id)
	if(is.null(atlas)) {
		stop("cannot find atlas id ", atlas_id, " when callig annotate mc by sc projection")
	}

	fig_cmp_dir = gen_fig_cmp_dir(fig_cmp_dir, qmc_id, atlas@atlas_name)

	a_mc = scdb_mc(atlas@mc_id)
	q_mc = scdb_mc(qmc_id)
	if(is.null(q_mc)) {
		stop("cannot find query mc id ", q_mc, " in scdb")
	}
	
	a_gset = scdb_gset(atlas@gset_id)

	a_gnms = names(a_gset@gene_set)
	a_gnms = intersect(a_gnms, rownames(a_mc@e_gc))

	q_egc = log2(q_mc@e_gc+1e-5)
	if(!is.null(q_gset_id)) {
		q_gset = scdb_gset(q_gset_id)
		q_gnms = names(q_gset@gene_set)
		q_egc = q_egc[q_gnms,]
	} 
	if(atlas@gene_naming_type != qmat_naming_type) {
		gene_name_src_targ = c(qmat_naming_type, atlas@gene_naming_type)
		q_egc = translate_mat_gnames(gene_name_src_targ, q_egc)
	}

	a_gnms = intersect(a_gnms, rownames(q_egc))
	
	q_egc = q_egc[a_gnms,]
	a_egc = log2(a_mc@e_gc[a_gnms,]+1e-5)

	a_self = tgs_cor_knn(a_egc, y=a_egc, knn=2)
	a_self = a_self[a_self$rank==2,]
	colnames(a_self) = c("a_mcid", "a_best", "a_top_cor", "rank")
	q_proj = tgs_cor_knn(q_egc, a_egc, knn=1)
	colnames(q_proj) = c("q_mcid", "a_mcid", "q_top_cor", "rank")
	q_proj = q_proj %>% left_join(a_self[,c("a_mcid","a_top_cor")])
	q_proj$dlt = q_proj$a_top_cor - q_proj$q_top_cor

	q_mc@colors[q_proj$q_mcid] = ifelse(q_proj$dlt < T_cor_gap, 
										a_mc@colors[q_proj$a_mcid], nd_color)

	if(!is.null(new_qmc_id)) {
		scdb_add_mc(new_qmc_id, q_mc)
	}

	png(sprintf("%s/mcmc_query_cor.png", fig_cmp_dir),w=600,h=600)
	plot(q_proj$q_top_cor, q_proj$a_top_cor, pch=19, 
				col=a_mc@colors[q_proj$a_mcid], 
				cex=0.6); 
	abline(a=0,b=1); 
	abline(a=T_cor_gap, b=1)
	dev.off()

	return(q_proj)
}

#' Annotate query metacell with atlas by comparing query cells to atlas MCs
#'
#' This will take each cell in the query MC and find its to correlated
#' metacell in the reference, then generating figures showing detailed comparison of how each metacell in the query is distributed in the atlas, and how the pool of all cells in the query MC compare to the pool of their best match projections
#'
#' @param atlas_id id of atlas object in scdb
#' @param qmat_id id of metacell object ina scdb
#' @param qmc_id id of metacell object ina scdb
#' @param qmat_naming_type naming scheme of query matrix/mc
#' @param fig_cmp_dir name of directory to put figures per MC 
#' @param q_gset_id query gene set id object (optional) - to restrict features used in comaprison to the intersection of atlas and query gsets
#' @param md_field metadata field too use as additional factor for plotting subsets
#' @param recolor_mc_id if this is specified,  the atlas colors will be projected on the query MCCs and updated to the scdb object  named recolor_mc
#'	@param plot_all_mc set this to T if you want a plot per metacell to show comparison of query pooled umi's and projected pooled umis.
#'
#' @export
#'

mcatlas_annotate_mc_by_sc2mc_projection = function(atlas_id,
				qmat_id, qmc_id, qmat_naming_type,
				fig_cmp_dir = NULL,
				q_gset_id = NULL,
				new_qmc_id = NULL,
				atlas_mc2d_id = NULL,
				plot_all_mcs = F,
				md_field=NULL,
				max_entropy=2,
				burn_cor=0.6)
{
	atlas = scdb_mcatlas(atlas_id)
	if(is.null(atlas_id)) {
		stop("atlas id ", atlas_id, " is missing")
	}
	a_mc = scdb_mc(atlas@mc_id)
	if(is.null(a_mc)) {
		stop("mc id ", atlas@mc_id, " is missing")
	}
	q_mc = scdb_mc(qmc_id)
	if(is.null(q_mc)) {
		stop("mc id ", mc_id, " is missing")
	}
	q_mat = scdb_mat(qmat_id)
	if(is.null(q_mat)) {
		stop("mat id ", qmat_id, " is missing")
	}

	a_gset = scdb_gset(atlas@gset_id)
	if(is.null(a_gset)) {
		stop("gset id ", atlas@gset_id, " is missing")
	}
	a_mc2d = NULL
	if(!is.null(atlas_mc2d_id)) {
		a_mc2d = scdb_mc2d(atlas_mc2d_id)
		if(is.null(a_mc2d)) {
			stop("atlas mc2d  ", atlas_mc2d_id, " is missing")
		}
	}

	a_gnms = names(a_gset@gene_set)
	a_gnms = intersect(a_gnms, rownames(a_mc@e_gc))

	q_umis = q_mat@mat[, names(q_mc@mc)]
	q_egc = q_mc@e_gc

	if(!is.null(q_gset_id)) {
		q_gset = scdb_gset(q_gset_id)
		q_gnms = names(q_gset@gene_set)
		q_umis = q_umis[q_gnms,]
		q_egc = q_egc[q_gnms,]
	} 
	if(atlas@gene_naming_type != qmat_naming_type) {
		gene_name_src_targ = c(qmat_naming_type, atlas@gene_naming_type)
		q_umis = translate_mat_gnames(gene_name_src_targ, q_umis)
		q_egc = translate_mat_gnames(gene_name_src_targ, q_egc)
	}

	a_gnms = intersect(a_gnms, rownames(q_umis))
	a_gnms = intersect(a_gnms, rownames(q_egc))
	
	q_umis = q_umis[a_gnms,]
	q_egc = q_egc[a_gnms,]

	a_legc = log2(a_mc@e_gc[a_gnms,]+1e-5)
	a_egc = a_mc@e_gc[a_gnms,]

	hits = tgs_cor_knn(as.matrix(q_umis), a_legc, knn=1)
	colnames(hits) = c("qsc_id", "amc_id", "cor", "rank")

	best_ref = hits$amc_id
	names(best_ref) = hits$qsc_id

	best_ref_cor = hits$cor
	names(best_ref_cor) = hits$qsc_id

	c_nms = names(q_mc@mc)
	mc_proj_col_p = table(q_mc@mc[c_nms], a_mc@colors[best_ref])
	mc_proj_col_p = mc_proj_col_p/rowSums(mc_proj_col_p)
	query_entropy = apply(mc_proj_col_p, 1, 
					function(x) { p =x/sum(x); lp =log2(p+1e-6); return(sum(-p*lp)) })
	
	if(!is.null(new_qmc_id)) {
		proj_col= a_mc@colors[as.numeric(unlist(best_ref))]
		new_col=tapply(proj_col,
						   q_mc@mc,
							function(x) { names(which.max(table(x))) }
							)
		q_mc@colors = new_col
		q_mc@colors[query_entropy > max_entropy] = "gray"
		scdb_add_mc(new_qmc_id, q_mc)
	}

	fig_cmp_dir = gen_fig_cmp_dir(fig_cmp_dir, qmc_id, atlas@atlas_name)

#plotting 2D porjection of query on atlas map
	if(!is.null(a_mc2d)) {
		png(sprintf("%s/comp_2d.png", fig_cmp_dir), w=1200, h=1200)
		layout(matrix(c(1,2),nrow=2), h=c(1,4))
		par(mar=c(0,3,2,3))
		plot(best_ref_cor[order(q_mc@mc)], pch=19, 
						col=a_mc@colors[best_ref[order(q_mc@mc)]],cex=0.6)
		grid()
		par(mar=c(3,3,0,3))
		n = length(best_ref)
		xrange = 0.02*(max(a_mc2d@mc_x) - min(a_mc2d@mc_x))
		yrange = 0.02*(max(a_mc2d@mc_y) - min(a_mc2d@mc_y))
		a_x = a_mc2d@mc_x[best_ref]+rnorm(n, 0, xrange)
		a_y = a_mc2d@mc_y[best_ref]+rnorm(n, 0, yrange)
		xlim = c(min(a_mc2d@mc_x), max(a_mc2d@mc_x))
		ylim = c(min(a_mc2d@mc_y), max(a_mc2d@mc_y))
		plot(a_x, a_y, pch=19, col=a_mc@colors[best_ref], ylim=ylim, xlim=xlim)
	}

#plotting best ref mc color distribution per mc
	atlas_cols = intersect(atlas@type_color,colnames(mc_proj_col_p))

	png(sprintf("%s/sc2mc_query_color_dist.png", fig_cmp_dir),w=600,h=150+12*nrow(mc_proj_col_p))
	layout(matrix(c(1,2),nrow=1),widths=c(5,2))
	par(mar=c(2,3,2,0))
	barplot(t(mc_proj_col_p[,atlas_cols]), col=atlas_cols, horiz=T, las=2)
	par(mar=c(2,0,2,3))
	barplot(query_entropy, col=q_mc@colors, horiz=T, las=2)
	grid()
	dev.off()

#comparing mean umis from projected mcs to mean umis on query mcs	
	q_mc_on_ref = t(tgs_matrix_tapply(a_egc[,best_ref], q_mc@mc, mean))
	rownames(q_mc_on_ref) = a_gnms

	cmp_lfp_1 = log2(1e-6+q_mc_on_ref)
	cmp_lfp_1n = cmp_lfp_1 - rowMeans(cmp_lfp_1)
	cmp_lfp_2 = log2(1e-6+q_egc)
	cmp_lfp_2n = cmp_lfp_2 - rowMeans(cmp_lfp_2)

	cross = tgs_cor(cmp_lfp_1n, cmp_lfp_2n)
	
	n = ncol(cmp_lfp_1)

	png(sprintf("%s/sc2mc_query_ref_cmp.png", fig_cmp_dir),w=600,h=600)
	shades = colorRampPalette(c("black", "darkblue", "white", "darkred", "yellow"))(1000)
	layout(matrix(c(1,4,2,3),nrow=2), h=c(10,1), w=c(1,10))
	par(mar=c(0,3,2,0))
	query_ref_colors = table(q_mc@mc, a_mc@colors[best_ref])
	query_mc_top_color = colnames(query_ref_colors)[apply(query_ref_colors,1,which.max)]
	image(t(as.matrix(1:length(q_mc@colors),nrow=1)), 
									col=query_mc_top_color, yaxt='n', xaxt='n')
	mtext(1:n,at=seq(0,1,l=n),side=2, las=1)
	par(mar=c(0,0,2,2))
	image(pmin(pmax(cross, -burn_cor),burn_cor), col=shades,xaxt='n', 
									yaxt='n', zlim=c(-burn_cor, burn_cor))
	par(mar=c(3,0,0,2))
	image(as.matrix(1:length(q_mc@colors),nrow=1), 
									col=q_mc@colors, yaxt='n', xaxt='n')
	mtext(1:n,at=seq(0,1,l=n),side=1, las=1)
	dev.off()
}

tba = function()
{

	if(plot_all_mcs) {
	ref_glob_p = rowMeans(ref_abs_fp)
	for(mc_i in 1:ncol(query_mc@mc_fp)) {
		ref_lfp = sort((query_mc_on_ref[,mc_i]+1e-5)/(ref_glob_p+1e-5))
		ref_marks = gene_name_map[names(tail(ref_lfp,15))]
		fig_nm = sprintf("%s/%d.png", fig_cmp_dir, mc_i)
		png(fig_nm, w=800, h=1200)
		layout(matrix(c(1,2), nrow=2))
	
		mcell_plot_freq_compare(query_mc@e_gc[common_genes,mc_i], 
						query_mc_on_ref[,mc_i], 
						n_reg = 1e-5,
						top_genes=20,
						highlight_genes = ref_marks,
						fig_h=600, fig_w=800,
						lab1="query", lab2="reference", 
						main=sprintf("reference/query, compare mc %d", mc_i))
		par(mar=c(2,2,2,2))
		plot(a_mc2d@sc_x, a_mc2d@sc_y, cex=0.2, pch=19, col="gray")
		f = query_mc@mc==mc_i
		points(ref_x[f], ref_y[f], pch=19, col=a_mc@colors[best_ref[f]])
		dev.off()
		
	}
	}
	if(!is.null(md_field)) {
		md = mat@cell_metadata
		if(md_field %in% colnames(md)) {
			md_vs = md[names(query_mc@mc),md_field]
			for(v in unique(md_vs)) {
				fig_nm = sprintf("%s/md_%s.png", fig_cmp_dir, v)
				png(fig_nm, w=800,h=800)
				f = (md_vs == v)
				plot(a_mc2d@sc_x, a_mc2d@sc_y, cex=0.2, pch=19, col="gray")
				points(ref_x[f], ref_y[f], pch=19, col=a_mc@colors[best_ref[f]])
				dev.off()
			}
		} else {
			message("MD field ", md_field, " not found")
		}
	}
	return(query_mc_on_ref)
	#plot compare_bulk
}


translate_mat_gnames = function(gene_names_src_targ, mat)
{
	gnames_df = scdb_gene_names_xref()
	key_s = gene_names_src_targ[1]
	key_t = gene_names_src_targ[2]
	f = !is.na(gnames_df[,key_t]) & !duplicated(gnames_df[,key_t])
	f = f & !is.na(gnames_df[,key_s]) & !duplicated(gnames_df[,key_s])
	src_2_targ= gnames_df[f,key_t]
	names(src_2_targ) = gnames_df[f, key_s]

	src_gnames = rownames(mat)
	f = src_gnames %in% names(src_2_targ) & 
						!duplicated(src_2_targ[src_gnames])

	if(sum(f) < 2) {
		stop("cannot translate ", length(gnames), " genes from ", key_s, " to ", key_t)
	}
	mat = mat[f,]
	rownames(mat) = src_2_targ[rownames(mat)]
		
	return(mat)
}

gen_fig_cmp_dir = function(fig_cmp_dir, qmc_id, atlas_id)
{
	if(is.null(fig_cmp_dir)) {
		fig_cmp_dir = paste(qmc_id, atlas_id, sep=".")
		fig_cmp_dir = paste(.scfigs_base, fig_cmp_dir,sep="/")
	}
	if(!dir.exists(fig_cmp_dir)) {
			dir.create(fig_cmp_dir)
	}
	return(fig_cmp_dir)
}
