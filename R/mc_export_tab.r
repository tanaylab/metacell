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
mcell_mc_export_tab = function(mc_id, gstat_id = NULL,
											T_gene_tot = 50,
											T_fold = 1,
											metadata_fields = c("batch_set_id"))
{
	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
		stop("MC-ERR non existing mc_id ", mc_id, " when trying to export fp table")
	}

	fp_max = apply(mc@mc_fp, 1, max)
	if(!is.null(gstat_id)) {
		stop("MC-ERR non existing gstat id ", gstat_id, " when trying to export fp table")
	}
	fp_tot = gstat[rowmames(mc@mc_fp),"tot"]
	
	# add #cells per clust and mean #umis (~ cell size)
	out_df = rbind(tapply(colSums(sc_cl@scmat@mat), sc_cl@clusts, mean), table(sc_cl@clusts))
	rownames(out_df) = c('mean_umis', 'n_cells')
	
	# add required breakdown to features
	if (!is.null(metadata_fields)) {
		for (s in metadata_fields) {
			out_df = rbind(table(sc_cl@scmat@cell_metadata[names(sc_cl@clusts), s], sc_cl@clusts), out_df)
		}
	}
	
	# actual clust_fp
	out_df = rbind(out_df, round(log2(sc_cl@clust_fp[f,]), 2))
	
	write.table(out_df, paste0(outdir, "/", tab_clust_fp_fn),	sep = "\t",quote = F)

}
