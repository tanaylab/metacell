#' colorize metacell using a set of prefered markers and their colors
#'
#' @param new_mc_id output metacell id in scdb
#' @param mc_id input metacell id in scdb (default: new_mc_id)
#' @param marker_color a data frame with fields gene, group, color, priority, thresh
#' @param override if this is true, all colors are going to be set to white unless some marker match is found
#'
#' @export
mc_colorize = function(new_mc_id, mc_id, marker_colors = NULL, override = T)
{
	sequential_coloring = get_param("mcp_colorize_by_seq_priority")
	
	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
		stop("MC-ERR metacell object is not avaialble in scdb, id = ", mc_id)
	}
	if(class(marker_colors)[1] != "data.frame"
	| length(intersect(c("gene","group", "color","priority","T_fold"),
								colnames(marker_colors))) != 5) {
		stop("MC-ERR marker colors parameter must be a data frame with fields gene, group, color, priority, T_fold")
	}
	marker_colors$gene = as.character(marker_colors$gene)
	marker_colors$color= as.character(marker_colors$color)
	rownames(marker_colors) = marker_colors$gene

	if(override) {
	  mc@colors = rep("white", ncol(mc@mc_fp))
	}
	good_marks = intersect(rownames(marker_colors),
						rownames(mc@mc_fp))
	if(length(good_marks) == 0) {
		message("no color markers are found")
		return
	}
	marker_colors = marker_colors[good_marks, ]

	mc@color_key = as.data.frame(marker_colors)

	cl_colors = rep(NA, ncol(mc@mc_fp))

	if (sequential_coloring) {
		for (p in sort(unique(marker_colors$priority))) {
			curr_marker_colors = marker_colors[marker_colors$priority == p, ]
			marker_fold = mc@mc_fp[curr_marker_colors$gene,]
			marker_fold = ifelse(marker_fold > curr_marker_colors$T_fold, marker_fold, NA)
			
			if (nrow(curr_marker_colors) == 1) {
				passed = is.na(cl_colors) & !is.na(marker_fold)
				hit = rep(1, sum(passed))
			}
			else {
				passed = is.na(cl_colors) & colSums(!is.na(marker_fold)) > 0
				hit = apply(marker_fold[, passed], 2, which.max)
			}
			
			cl_colors[passed] = curr_marker_colors[hit, 'color']
		}
	} else {
		marker_colors = marker_colors[order(marker_colors$priority),]
		marker_fold = mc@mc_fp[marker_colors$gene,]
		marker_fold = ifelse(marker_fold > marker_colors$T_fold, log2(marker_fold), NA)
		marker_fold = marker_fold * marker_colors$priority

		if(length(good_marks) > 1) {
			nonz = colSums(!is.na(marker_fold)) > 0
			hit = apply(marker_fold[, nonz], 2, which.max)
		} else {
			nonz = marker_fold > 0
			hit = rep(1, sum(nonz))
		}

		cl_colors[nonz] = marker_colors[hit, "color"]
	}
	
	if(!override) {
		cl_colors[is.na(cl_colors)] = mc@colors[is.na(cl_colors)]
	}
	mc@colors = cl_colors
	scdb_add_mc(new_mc_id, mc)
}

#' colorize metacell using an ugly default color spectrum, or a user supplied one
#'
#' @param mc_id metacell id in scdb
mc_colorize_default = function(mc_id, spectrum = NULL)
{
	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
		stop("MC-ERR metacell object is not avaialble in scdb, id = ", mc_id)
	}

	if(is.null(spectrum)) {
		spectrum = colorRampPalette(c("white", "lightgray", "darkgray", "burlywood1", "chocolate4","orange", "red", "purple", "blue", "cyan"))
	}
	mc@colors = spectrum(max(mc@mc))
	scdb_add_mc(mc_id, mc)
}
