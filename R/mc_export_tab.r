#' Output genes cell modules footprint matrix with metadata on the cell modules
#'
#' @param mc_id meta cell object id
#' @param gstat_id gstat object id
#' @param T_gene_tot minimal number of umis to include gene (to reduce table size)
#' @param T_fold threshold for maximal fold change over metacells
#' @param metadata_fields  metadata field names to be used as factors and breakdown over mcs
#'
#' @return Nothing. But save a table in the figure directory that can be loaded to excel
#'
#' @export
#'
mcell_mc_export_tab = function(mc_id, gstat_id = NULL, mat_id = NULL,
											T_gene_tot = 50,
											T_fold = 1, # should we also export genes with depletion (i.e. below certain threshold)?
											metadata_fields = c("batch_set_id"))
{
	mc = scdb_mc(mc_id)
	gstat = scdb_gstat(gstat_id)
	scmat = scdb_mat(mat_id)
	if(is.null(mc)) {
		stop("MC-ERR non existing mc_id ", mc_id, " when trying to export fp table")
	}

	if(is.null(gstat)) {
		stop("MC-ERR non existing gstat id ", gstat_id, " when trying to export fp table")
	}

	if(is.null(scmat)) {
		stop("MC-ERR non existing mat id ", mat_id, " when trying to export fp table")
	}

	# add #cells per clust, mean #umis (~ cell size) and assigned group (if exist)
	if (nrow(mc@color_key) > 0) {
		col2group = as.character(mc@color_key$group)
		names(col2group) = as.character(mc@color_key$color)

		groups = col2group[mc@colors]
	}
	else {
		groups = rep(NA, max(mc@mc))
	}

	out_df = rbind(tapply(colSums(scmat@mat[, names(mc@mc)]), mc@mc, mean), table(mc@mc), groups, seq_along(groups))
	out_df = cbind(rep("", nrow(out_df)), out_df)
	rownames(out_df) = c('mean_umis', 'n_cells', 'group', 'mc_id')

	# add required breakdown to features
	if (!is.null(metadata_fields)) {
		for (s in metadata_fields) {
		    new_df = table(scmat@cell_metadata[names(mc@mc), s], mc@mc)
		    new_df = cbind(rep(s, nrow(new_df)), new_df)
			out_df = rbind(new_df, out_df)
		}
	}

	fp_max = apply(mc@mc_fp, 1, max)
	fp_tot = gstat[intersect(rownames(mc@mc_fp), rownames(gstat)), "tot"]

	# genes to export
	f = fp_max > T_fold & fp_tot > T_gene_tot

	# actual clust_fp
	out_df = rbind(out_df, cbind(rep("", sum(f)), round(log2(mc@mc_fp[f,]), 2)))

	out_df = cbind(out_df[, 1], rownames(out_df), out_df[, -1])
	tab_clust_fp_fn = sprintf("%s/%s.log2_mc_fp.txt", .scfigs_base, mc_id)
	write.table(out_df, tab_clust_fp_fn, sep = "\t", quote = F, row.names = F, col.names = F)
}
