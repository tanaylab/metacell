
#' Create a gene name xref for heuristc map of mars/10x
#'
#'
#'

mcell_add_gene_names_xref_mars10x = function(mars_mat_id, tenx_mat_id, xref_id="DB")
{
	mat1 = scdb_mat(mars_mat_id) 
	if(is.null(mat1)) {
		stop("no mars mat ", mars_mat_id)
	}
	mat2 = scdb_mat(tenx_mat_id) 
	if(is.null(mat2)) {
		stop("no tenx mat ", tenx_mat_id)
	}

	mars_m_nms = rownames(mat1@mat)

	mars_nms = strsplit(mars_m_nms, ";")
	mars_n = unlist(lapply(mars_nms, length))
	mars_v = rep(mars_m_nms,times=mars_n)
	names(mars_v) = unlist(mars_nms)

	tenx_m_nms = rownames(mat2@mat)

	df1 = data.frame(key = names(mars_v), mars = mars_v)
	df2 = data.frame(key = tenx_m_nms, tenx = tenx_m_nms)
	xref_df = df1 %>% full_join(df2)

	scdb_add_gene_names_xref(xref_id, xref_df)
}

#' Generate mapping of 10x to mars names
#'
#'  Not more than finding which concatenated names (";" delimited) are related to 10x gene names. Should be replaced by something more systematic that will happen during import
#'
#' @param mars_mc_id metacell id of a mars-seq dataset
#' @param tenx_mc_id metacell id of a 10x dataset
#'
#' @export
gen_10x_mars_gene_match = function(mars_mc_id, tenx_mc_id)
{
	mc1 = scdb_mc(mars_mc_id) 
	if(is.null(mc1)) {
		stop("not mars mc ", mars_mc_id)
	}
	mc2 = scdb_mc(tenx_mc_id) 
	if(is.null(mc2)) {
		stop("not tenx mc ", tenx_mc_id)
	}

	mars_m_nms = rownames(mc1@e_gc)

	mars_nms = strsplit(mars_m_nms, ";")
	mars_n = unlist(lapply(mars_nms, length))
	mars_v = rep(mars_m_nms,times=mars_n)
	names(mars_v) = unlist(mars_nms)

	nm10x = rownames(mc2@e_gc)
	ten2mars = mars_v[nm10x]

	f = !duplicated(ten2mars)
	mars2ten = nm10x[f]
	names(mars2ten) = ten2mars[f]

	return(ten2mars)
}
