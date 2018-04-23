#' build a metacell cover from a co-clust object using a simple hclust approach
#'
#' @param mc_id Id of new metacell object
#' @param coc_id cocluster object to use
#' @param mat_id mat object to use when building the mc object
#' @param K number of clusters to generate from the coclust hclust. Look at the coclust plot to detemrine this. But make sure you keep the clusters at reasonable size (e.g. <100) for use as metacells.
#' @param force set this to true to overide size limits on hclust
#'
mcell_mc_from_coclust_hc = function(mc_id, coc_id, mat_id, K, force=F)
{
	tgs_coclust_hc_type = get_param("scm_coclust_hc_type")

	coc = scdb_coclust(coc_id)
	if(is.null(coc)) {
		stop("MC-ERR: coclust ", coc_id , " is missing when running mc from coclust")
	}
	mat = scdb_mat(mat_id)
	if(is.null(mat)) {
		stop("MC-ERR: mat id ", mat_id, " is missing when running add_mc_from_graph")
	}
	if(ncol(mat@mat) > 20000 & !force) {
		stop("no support for hclust on coclust with more than 20K nodes, use force=T")
	}
	message("running cocluster hclust now")

	m_coc = as.matrix(sparseMatrix(coc@coclust$node1, coc@coclust$node2, x=coc@coclust$cnt))
	m_samp = (coc@n_samp/mean(coc@n_samp)) %*% t(coc@n_samp)
	m_coc = m_coc/m_samp
	rownames(m_coc) = colnames(mat@mat)
	colnames(m_coc) = colnames(mat@mat)
	outliers = colnames(mat@mat)[rowSums(m_coc)==0]
	m_coc = m_coc[rowSums(m_coc) > 0,rowSums(m_coc)>0]
	hc = hclust(as.dist(1-tgs_cor(m_coc)), tgs_coclust_hc_type)
	clust = cutree(hc, K)
	scdb_add_mc(mc_id, tgMCCov(clust, outliers, mat))
	message("reordering metacells by hclust and most variable two markers")
	mcell_mc_reorder_hc(mc_id)
}

#' build a metacell cover from a big co-clust using louvain clustering and metacell coverage within clusters
#'
#' @param mc_id Id of new metacell object
#' @param coc_id cocluster object to use
#' @param mat_id mat object to use when building the mc object
#' @param max_mc_size maximum mc size (bigger clusters will be dissected)
#' @param max_clust_size maximum clust size. Bigger chunks will be clustered
#'
mcell_mc_from_coclust_louv_sub = function(mc_id, coc_id, mat_id,
		max_clust_size, max_mc_size, min_mc_size, T_weight = 1)
{

#	tgs_coclust_hc_type= get_param("scm_coclust_hc_type")

	coc = scdb_coclust(coc_id)
	if(is.null(coc)) {
		stop("MC-ERR: coclust ", coc_id , " is missing when running mc from coclust")
	}
	mat = scdb_mat(mat_id)
	if(is.null(mat)) {
		stop("MC-ERR: mat id ", mat_id, " is missing when running add_mc_from_graph")
	}

#hierarchically break coclust
	filt_coc = coc@coclust[coc@coclust$cnt > T_weight,]
	tot_deg = tabulate(c(filt_coc$node1, filt_coc$node2), nbins = ncol(mat@mat))
	outliers = colnames(mat@mat)[tot_deg <= 1]
	h_mc = list(root=which(tot_deg>1))
	h_coc = list(root=filt_coc)
	recent_zoom=TRUE
	while(recent_zoom) {
		new_h_mc = list()
		new_h_coc = list()
		recent_zoom = FALSE
		message("next h iteration, current sup_clst N = ", length(h_mc))
		for(cl_i in 1:length(h_mc)) {
			next_id = length(new_h_mc) + 1
			sz = length(h_mc[[cl_i]])
			if(sz > max_clust_size) {
				sub_clust = coc_dissect_louvain(h_coc[[cl_i]], h_mc[[cl_i]])
				message("open up large cluster on ", sz, " nodes, got ", length(sub_clust), " sub clusts, sizes ", paste(lapply(sub_clust, length),collapse=","))
				f = lapply(sub_clust, length) < min_mc_size
				if(sum(f) > 0) {
					outliers = c(outliers, colnames(mat@mat)[unlist(sub_clust[f])])
					message("removing ", sum(f), " small clusters into the outlier set, new total outliers ", length(outliers))
					sub_clust = sub_clust[!f]
				}
				if(length(sub_clust) > 0) {
					names(sub_clust) = next_id : (next_id + length(sub_clust) - 1)
					new_h_mc = append(new_h_mc, sub_clust)
					sub_coc = gen_sub_coclust(h_coc[[cl_i]],sub_clust)
					new_h_coc = append(new_h_coc, sub_coc)
				}

				#break using louvain
				if(length(sub_clust) > 1) {
					recent_zoom = TRUE
				} else {
					message("louvain returned just one cluster on N =", length(h_mc[[cl_i]]), "\n")
				}
			} else if(sz > max_mc_size) {
				sub_mc = coc_dissect_mc(h_coc[[cl_i]], h_mc[[cl_i]], min_mc_size)
				message("cover intermediate cluster on ", sz, " nodes, got ", length(sub_mc), " new mcs, size ", paste(lapply(sub_mc, length),collapse=","))
				if(length(sub_mc) == 1 & length(sub_mc[["0"]])>0) {
					message("renaming ", length(sub_mc[["0"]]), " outliers, tot edges was ", nrow(h_coc[[cl_i]]))
					names(sub_mc) = c(1)
				}
				if(length(sub_mc[["0"]])>0) {
					message("adding ", length(sub_mc[["0"]]), " outliers, tot edges was ", nrow(h_coc[[cl_i]]))
					outliers = c(outliers, colnames(mat@mat)[sub_mc[["0"]]])
					sub_mc[["0"]] = NULL
				}
				if(length(sub_mc) != 0) { #only if all of it is outliers
					names(sub_mc) = next_id : (next_id + length(sub_mc) - 1)
					new_h_mc = append(new_h_mc, sub_mc)
					sub_coc = gen_sub_coclust(h_coc[[cl_i]],sub_mc)
					new_h_coc = append(new_h_coc, sub_coc)
				} else {
					browser()
				}
			} else {
				message("retain mc ", cl_i, " as ", next_id, " N = ", length(h_mc[[cl_i]]))
				new_h_mc[[next_id]] = h_mc[[cl_i]]
				new_h_coc[[next_id]] = h_coc[[cl_i]]
				names(new_h_mc)[[next_id]] = next_id
				names(new_h_coc)[[next_id]] = next_id
			}
		}
	   h_mc = new_h_mc
		h_coc = new_h_coc
	}

	clust = rep(as.numeric(names(h_mc)), times=lapply(h_mc, length))
	names(clust) = colnames(mat@mat)[unlist(h_mc)]
	scdb_add_mc(mc_id, tgMCCov(clust, outliers, mat))
	message("reordering metacells by hclust and most variable two markers")
	mcell_mc_reorder_hc(mc_id)
}


coc_dissect_louvain = function(coclust, nodes)
{
	colnames(coclust) = c("node1", "node2", "weight")

	message("will build graph")
	g_cc = graph_from_data_frame(d=coclust, directed=F)
	message("will cluster w louvain")
	a = cluster_louvain(g_cc)

	clusts = split(nodes, a$membership)
	return(clusts)
}

coc_dissect_mc = function(coclust, nodes, min_mc_size)
{
	tgs_clust_cool = get_param("scm_tgs_clust_cool")
	tgs_clust_burn = get_param("scm_tgs_clust_burn_in")

	edges = coclust
	colnames(edges) = c("col1", "col2", "weight")
	edges$weight = edges$weight/max(edges$weight)

	edges = edges[edges$col1 != edges$col2,]
	redges = edges[,c("col1", "col2", "weight")]
	redges$col1= edges$col2
	redges$col2= edges$col1
	edges = rbind(edges,redges)

	node_clust = tgs_graph_cover(edges, min_mc_size,
					cooling = tgs_clust_cool, burn_in = tgs_clust_burn)
	rownames(node_clust) = node_clust$node
	clusts = split(nodes, node_clust$cluster[nodes])
	return(clusts)
}

gen_sub_coclust = function(coc, subc)
{
	message("dissect ", nrow(coc), " edges ")
	clust = rep(as.numeric(names(subc)), times=lapply(subc, length))
	clust_map = rep(NA, max(unlist(subc)))
	clust_map[unlist(subc)] = clust
	c1 = clust_map[coc$node1]
	c2 = clust_map[coc$node2]
	message("done building map")
	f = !is.na(c1) & !is.na(c2) & c1 == c2
	message("done getting factor")
	coc = coc[f,]
	intra_c = c1[f]
	message("done first subseting")
	names(intra_c) = NULL
	sub_coc = split(coc, intra_c)
	names(sub_coc) = names(subc)
	return(sub_coc)
}

#' build a metacell cover from co-clust data through filtering un-balanced edges
#' and running graph cover 
#'
#' @param mc_id Id of new metacell object
#' @param coc_id cocluster object to use
#' @param mat_id mat object to use when building the mc object
#' @param K - this will 
#' @param min_mc_size minimum mc size for graph cov
#' @param alpha the threshold for filtering edges by their coclust weight is alpha * (Kth highest coclust on either node1 or node2)
#'
mcell_mc_from_coclust_balanced = function(mc_id, coc_id, 
												mat_id, K, min_mc_size, alpha=2)
{
	
	tgs_clust_cool = get_param("scm_tgs_clust_cool")
	tgs_clust_burn = get_param("scm_tgs_clust_burn_in")

	coc = scdb_coclust(coc_id)
	if(is.null(coc)) {
		stop("MC-ERR: coclust ", coc_id , " is missing when running mc from coclust")
	}
	mat = scdb_mat(mat_id)
	if(is.null(mat)) {
		stop("MC-ERR: mat id ", mat_id, " is missing when running add_mc_from_graph")
	}

#hierarchically break coclust
	edges = coc@coclust
	deg_wgt = as.matrix(table(c(edges$node1, edges$node2), c(edges$cnt,edges$cnt)))
	deg_cum = t(apply(deg_wgt, 1, function(x) cumsum(rev(x))))
	thresh_Kr = rowSums(deg_cum > K)
	thresh_K = rep(NA, length(levels(edges$node1)))
	names(thresh_K) = levels(edges$node1)
	thresh_K[as.numeric(names(thresh_Kr))] = thresh_Kr

	filt_edges = thresh_K[edges$node1] < edges$cnt * alpha | 
							thresh_K[edges$node2] < edges$cnt * alpha

	message("filtered ", nrow(edges) - sum(filt_edges), " left with ", sum(filt_edges), " based on co-cluster imbalance")
	edges = edges[filt_edges,]

	colnames(edges) = c("col1", "col2", "weight")
	edges$weight = edges$weight/max(edges$weight)

	edges = edges[edges$col1 != edges$col2,]
	redges = edges[,c("col1", "col2", "weight")]
	redges$col1= edges$col2
	redges$col2= edges$col1
	edges = rbind(edges,redges)

	node_clust = tgs_graph_cover(edges, min_mc_size, 
					cooling = tgs_clust_cool, burn_in = tgs_clust_burn)

	f_outlier = (node_clust$cluster == 0)

	outliers = colnames(mat@mat)[node_clust$node[f_outlier]]
	mc = as.integer(as.factor(node_clust$cluster[!f_outlier]))
	names(mc) = colnames(mat@mat)[!f_outlier]
	message("building metacell object, #mc ", max(mc))
	cell_names = colnames(mat@mat)
	scdb_add_mc(mc_id, tgMCCov(mc, outliers, mat))
	message("reordering metacells by hclust and most variable two markers")
	mcell_mc_reorder_hc(mc_id)
}
