#' generate a batch stat table for a given single cell umi matrix
#' @param scmat - the matrix to work on
#' @return batch stats table
#'
#' @export

mcell_batch_stats =  function(scmat)
{
	md = scmat@cell_metadata
	umis = scmat@mat
	# report batch stats
	batch_stats = data.frame()
	batch_nms = sort(unique(md$amp_batch_id))
	for(b in batch_nms) {
		batch_stats = rbind(batch_stats, mcell_calc_one_batch_stats(scmat, b))
	}
	rownames(batch_stats) = batch_nms

	return(batch_stats)
}
#	write.table(batch_stats, sprintf("%s/batch_stats.txt", get_param("outdir")), quote=F, sep="\t")
#	.plot_batch_stats(batch_stats, batch_meta , sprintf("%s/batch_stats.png", get_param("outdir")))


#' calc batch stats - essentially umi distribution
#'
#' @param scmat - a matrix object
#' @param batch_id - the batch name to summarize
#'
#' @return one row data frame for merging
#' @export
#'
mcell_calc_one_batch_stats = function(scmat, batch_id)
{
	md = scmat@cell_metadata
	mat = scmat@mat

	f_batch = md$amp_batch_id == batch_id

	n_c = colSums(mat[,f_batch])

	data.frame(
		total = sum(n_c),
		total_ercc = sum(md$spike_count[f_batch]),
		n_lt50 = sum(n_c < 50),
		n_50_100 = sum(n_c >=  50 & n_c < 100),
		n_100_200 = sum(n_c >=  100 & n_c < 200),
		n_200_400 = sum(n_c >=  200 & n_c < 400),
		n_400_800 = sum(n_c >=  400 & n_c < 800),
		n_800_1600 = sum(n_c >=  800 & n_c < 1600),
		n_gt1600 = sum(n_c >= 1600),
		n_valid = length(n_c),
		valid_q25 = quantile(n_c, 0.25),
		valid_q50 = quantile(n_c, 0.5),
		valid_q75 = quantile(n_c, 0.75)
		)
}

#' plot batches stats
#'
#' @param batches_stats collected by .calc_batch_stats
#' @param fig_nm output file name
#' @param cells_breakdown_field color #cells per batch by this field
#' @param cell_col_dict color dictionary for cells_breakdown_field values
#'
#' @export
#'
mcell_plot_batch_stats = function(mat_id, fig_nm = NULL)
{
	scmat = scdb_mat(mat_id)
	if(is.null(fig_nm)) {
		fig_nm = scfigs_fn(mat_id, "batchstat")
	}

	cells_breakdown_field=get_param("scm_cells_breakdown_field")
	cell_col_dict=get_param("scm_cell_col_dict")

	batch_stats = mcell_batch_stats(scmat)
	md = scmat@cell_metadata
	batch_meta = unique(md[, grep("spike_count", colnames(md), invert=T)])

	rownames(batch_meta) = batch_meta$amp_batch_id

	batch_stats = batch_stats[order(gsub("AB", "", rownames(batch_stats)), decreasing = T), ]
	batch_meta = batch_meta[rownames(batch_stats), ]

	b_tot = batch_stats$total
	names(b_tot) = rownames(batch_stats)

	main_cex = 2
	bp_cex = 1

	.plot_start(fig_nm, w=960, h=max(800, 100 + 10 * length(b_tot)))
	layout(matrix(1:4, 1, 4), w=c(2,1,5,1))

	type_cols = 'blue'
	if (!is.null(cell_col_dict) & !is.null(cells_breakdown_field)) {
		type_cols = unlist(cell_col_dict[batch_meta[, cells_breakdown_field]])
	}

	par(mar=c(4,5,2,0))
	barplot(b_tot, horiz=T, las=2, cex.axis=bp_cex, cex.names=bp_cex, border=NA, col=type_cols)
	title(main='total UMIs', cex.main=main_cex)
	if (!is.null(cell_col_dict) & !is.null(cells_breakdown_field)) {
		legend("topleft", legend=names(type_cols), bty='n', border=NA, fill=unlist(cell_col_dict), cex=bp_cex)
	}
	box(lwd=0.5)

	par(mar=c(4,1,2,0))
	barplot(batch_stats$total_ercc/batch_stats$total, cex.axis=bp_cex, horiz=T, border=NA, col='black')
	title(main='%ERCC', cex.main=main_cex)
	box(lwd=0.5)

	bd_cols = RColorBrewer::brewer.pal(n=7, 'Blues')
	barplot(t(batch_stats[,3:9]), horiz=T, border=NA, cex.axis=bp_cex, col=bd_cols, yaxt='n')
	title(main='#cells by #UMIs', cex.main=main_cex)
	legend("topleft", legend=gsub("_", "-", gsub("n_", "", colnames(batch_stats)[3:9])), bty='n', fill=bd_cols, cex=1.2 * bp_cex, pt.cex=1.2 * bp_cex, ncol=length(bd_cols), border=NA, x.intersp=0.8)
	box(lwd=0.5)

	par(mar=c(4,1,2,1))
	barplot(batch_stats$n_valid, horiz=T, border=NA, cex.axis=bp_cex, col='black')
	title(main='#valid', cex.main=main_cex)
	box(lwd=0.5)

	dev.off()
}

#' Plot histogram of total number of umis per cell in the umis matrix
#'
#' @param mat_id id of the mat to use
#' @param min_umis_cutoff If not null, will mark this cutoff by a line. If NA, auto calculate by finding the first bin in the histogram with more cells then its predecessor.
#' @param bin_for_cutoff bin size for histogram used to auto find min_umis_cutoff
#'
#' @return min_umis_cutoff (either the given, the calculated, or null)
#'
#' @export
#'
mcell_plot_umis_per_cell = function(mat_id, min_umis_cutoff = NA, bin_for_cutoff = 50)
{
  mat = scdb_mat(mat_id)

  if (is.null(mat)) {
    stop(sprintf("MCERR: mat with id %s not found", mat_id))
  }

  uc = Matrix::colSums(mat@mat)
  .plot_start(scfigs_fn(mat_id, "total_umi_distr"), w=450, h=450)

  if (is.na(min_umis_cutoff)) {
    h = hist(uc, breaks=seq(0, max(uc) + bin_for_cutoff, by=bin_for_cutoff), plot=F)
    min_umis_cutoff = min(which(diff(h$counts) > 0)) * bin_for_cutoff
  }
  h = hist(log2(uc), 200, xlab='total cell UMIs (log2)', main=mat_id)

  if (!is.null(min_umis_cutoff)) {
    abline(v=log2(min_umis_cutoff), col='red', lty=2)
    text(log2(min_umis_cutoff), max(h$count)/2, min_umis_cutoff, pos=2, col='red')
  }

  dev.off()

  return(min_umis_cutoff)
}
#
