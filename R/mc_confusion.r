
#' Compute confusion matrix on metacells
#'
#' @param mc_id meta cell object
#' @param graph_id graph object
#' @param K Knn parameter for filtering graph edges
#'
#' @export

mcell_mc_confusion_mat = function(mc_id, graph_id, K, ignore_mismatch = F)
{
	graph = scdb_cgraph(graph_id)
	if(is.null(graph)) {
		stop("MC-ERR: cell graph id ", graph_id, " is missing when computing mc confusion")
	}
	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
		stop("MC-ERR: mc id ", mc_id, " is missing when computing mc confusion")
	}
	if(!ignore_mismatch & !mcell_mc_match_graph(mc_id, graph_id)) {
		stop("MC-ERR: graph is missing some of the cells in the mc cover when computing confusion, halting")
	}

#filter top K edges
	deg = nrow(graph@edges)/length(graph@nodes)
	T_w = 1-(K+1)/deg
	edges = graph@edges
	f = edges$w > T_w
#factor of edges$mc1 and factor of mc may be different
	mc_map = mc@mc[as.character(levels(edges$mc1))] 
	mc1 = mc_map[edges$mc1[f]]
	mc2 = mc_map[edges$mc2[f]]
	N = max(mc1, mc2, na.rm=T)

	confu = matrix(tabulate((mc1-1)*N+mc2, nbins=N*N), nrow=N)
	rownames(confu) = 1:N
	colnames(confu) = 1:N
	return(confu)
}

#' Compute cell homogeneity - the fraction of intra mc edges per cell
#'
#' @param mc_id meta cell object
#' @param graph_id graph object
#'
#' @export

mcell_mc_cell_homogeneity = function(mc_id, graph_id, ignore_mis=F)
{
	graph = scdb_cgraph(graph_id)
	if(is.null(graph)) {
		stop("MC-ERR: cell graph id ", graph_id, " is missing when computing mc confusion")
	}
	mc = scdb_mc(mc_id)
	if(is.null(mc)) {
		stop("MC-ERR: mc id ", mc_id, " is missing when computing mc confusion")
	}
	if(!ignore_mis & !mcell_mc_match_graph(mc_id, graph_id)) {
		stop("MC-ERR: graph is missing some of the cells in the mc cover when computing confusion, halting")
	}
#filter top K edges
	edges = graph@edges
#factor of edges$mc1 and factor of mc may be different
	mc_map = mc@mc[as.character(levels(edges$mc1))] 
	intra_mc = mc_map[edges$mc1] == mc_map[edges$mc2]
	f = !is.na(intra_mc)

	N = max(as.integer(edges$mc1),na.rm=T)
	intra_stat = tabulate(as.integer(edges$mc1[f])*2 - ifelse(intra_mc[f],0,1), nbins = N*2)

	homogeneity = t(matrix(intra_stat, nrow=2))
	rownames(homogeneity) = levels(edges$mc1)
	homogeneity = homogeneity[intersect(levels(edges$mc1), names(mc@mc)),]
	colnames(homogeneity) = c("out", "in")
	return(homogeneity)
}
