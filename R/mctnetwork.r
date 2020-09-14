#' temporal netwrok over metacells
#'
#' Splitting metacells over a discrete time axis, defining manifold connections and estimated flows over them
#'
#' @slot mc_id id of the metacell object we represent as a network
#' @slot mgraph_id id of the mgraph object defining manifold structure
#' @slot times_nms names of the time points (Default 1:T)
#' @slot mc_t distribution of metacells (rows) over time points (cols)
#' @slot mc_mgraph a data frame defining triplets mc1, mc2, distance.
#' @slot mc_manif_p probability of moving between mcs - matrix
#' @slot mc_manif_cost  probability of moving between mcs
#' @slot network - a data frame defining the network structure
#' @slot edge_flows - flow of mincost solution
#' @slot mc_t_infer - inferred fraction ofcells per mc and time
#'
#' @export tgMCTNetwork
#' @exportClass tgMCTNetwork
tgMCTNetwork <- setClass(
   "tgMCTNetwork",
	slots = c(
	  mc_id = "character",
	  mgraph_id = "character",
	  cell_time = "vector",
	  mc_t = "matrix",
	  mc_mgraph = "data.frame",
	  mc_cost_mat = "matrix",
	  mc_manif_cost = "data.frame",
	  network = "data.frame",
	  edge_flows = "vector",
	  mc_t_infer = "matrix",
	  mc_forward = "list",
	  mc_backward = "list")
)

#' Construct a meta cell time network
#'
#'
#'
#' @param mc_id metacell object id 
#' @param mgraph_id metacell graph object id 
#' @param cell_time vector assigning a time (number) to each cell
#' @param mc_prolif factor per mc
#' @export

setMethod(
  "initialize",
  signature = "tgMCTNetwork",
  definition =
    function(.Object, mc_id, mgraph_id, cell_time) {
		.Object@mc_id = mc_id
		.Object@mgraph_id = mgraph_id
		.Object@cell_time = cell_time
		mc = scdb_mc(mc_id)
		if(is.null(mc)) {
			stop("MC-ERR unkown mc_id ", mc_id, " when building network")
		}
		mc_t = table(mc@mc, cell_time[names(mc@mc)])
		mc_t = matrix(mc_t, ncol=ncol(mc_t))
		mc_t = t(t(mc_t)/colSums(mc_t))

		.Object@mc_t = matrix(mc_t, ncol=ncol(mc_t))

		mgraph = scdb_mgraph(mgraph_id)
		if(is.null(mgraph)) {
			stop("MC-ERR unkown manifold id ", mgraph_id, " when building network")
		} 
		if(mgraph@mc_id != mc_id) {
			stop("MC-ERR mismatch of mc_id ", mc_id, " and mgraph ", mgraph_id, " when building network")
		}
		.Object@mc_mgraph = mgraph@mgraph

      return(.Object)
    }
)

#' Generate a new network in scdb
#'
#' This constructs a meta cell time network object from an mc and mgraph objects and assignments of cells to graphs
#'
#' @param net_id id of scdb network object ot be added
#' @param mc_id metacell object
#' @param mgraph_id metacell manfold graph object
#' @param cell_time assigning of time (discrete) to cells
#' @export
mcell_new_mctnetwork = function(net_id, mc_id, mgraph_id, cell_time)
{
	scdb_add_mctnetwork(net_id, tgMCTNetwork(mc_id, mgraph_id, cell_time))
}

#' Compute and update manifold costs for a mc time network
#'
#'
#' @param mct mct network object
#' @param t_exp time parameter for the markov exponential
#' @param T_cost threshold on the cost for retan
#'
#' @return an updated mct object
#' @export

mctnetwork_comp_manifold_costs = function(mct, t_exp = 1, T_cost = 1000)
{
  mgraph = mct@mc_mgraph	
  df = mgraph[,c(2,1,3)]
  colnames(df) = colnames(mgraph)
  
  mgraph2 = rbind(mgraph, df)
  mgraph2 = unique(mgraph2[,c(1,2,3)])
  mgraph2$dist_inv = 1/mgraph2$dist
  
  adj_mat = sparseMatrix(i = mgraph2$mc1,j = mgraph2$mc2,x = mgraph2$dist_inv)
  adj_mat = as.matrix(adj_mat)
  
  diag(adj_mat) = 0
  
  row_max = apply(adj_mat,1,max)
  median_rate = median(row_max)
  
  adj_mat = adj_mat/median_rate
  
  row_sums = rowSums(adj_mat)
  diag(adj_mat) = -row_sums
  
  trans_mat = expm(t_exp*adj_mat)
  
  trans_mat = as.matrix(trans_mat)
  trans_mat = trans_mat/apply(trans_mat,1,max)

#  diag(trans_mat) = rowMaxs(trans_mat)

  cost_mat = round(10/trans_mat)
  cost_mat = as.matrix(cost_mat)
  rownames(cost_mat) = 1:nrow(cost_mat)
  colnames(cost_mat) = 1:ncol(cost_mat)
  
  manifold = as.data.frame.table(cost_mat)
  colnames(manifold) = c("mc1","mc2","cost")
  manifold = manifold[is.finite(manifold$cost),]
  manifold = manifold[manifold$cost < T_cost,]
  
  colnames(manifold) = c("mc1","mc2","cost")

	mct@mc_cost_mat = cost_mat
	mct@mc_manif_cost = manifold
	return(mct)
}

#' Compute and update network strcuture from mc, manifodld costs
#'
#' @param mct mct network object
#' @param mc_leak the fraciton of frequency loss per mc per time step. All zero if you don't know what to expect.
#'
#' @return an updated mct object
#' @export

mctnetwork_gen_network = function(mct, mc_leak,
								capacity_var_factor = NULL,
								off_capacity_cost1 = 1,
								off_capacity_cost2 = 1000,
								k_norm_ext_cost = 2,
								k_ext_norm_cost = 2,
								k_ext_ext_cost = 100,
								T_cost = 1000) 
{
	k_inf = 100000

	manifold = mct@mc_manif_cost
	mc_t_freq = mct@mc_t
	max_t = ncol(mc_t_freq)
	gl = build_growth_mats(mct, mc_leak)
	growth_leak = gl$leak
	growth_compensate = gl$compensate

	if(is.null(capacity_var_factor)) {
		capacity_var_factor = rep(0.25, nrow(mc_t_freq))
	}

	edges = data.frame(from = c(), to = c(), 
						ID = c(), capacity = c(), cost = c(), 
						time1 = c(), time2 = c(), min_capacity = c(), 
						type1 = c(), type2 = c())
#type are: norm, ext, growth, src ,sink

	n_mc = nrow(mc_t_freq)

#generate lists of nodes per network layer

	source_id = 1

	normal_nodes_t_back = list()
	extend_nodes_t_back = list()	
	normal_nodes_t_front = list()	
	extend_nodes_t_front = list()	
	growth_nodes_t = list()
	next_n_id = 2
	for(t in 1:max_t) {
		normal_nodes_t_back[[t]] = next_n_id:(next_n_id + n_mc - 1)
		next_n_id = next_n_id + n_mc
		extend_nodes_t_back[[t]] = next_n_id:(next_n_id + n_mc - 1)
		next_n_id = next_n_id + n_mc
		normal_nodes_t_front[[t]] = next_n_id:(next_n_id + n_mc - 1)
		next_n_id = next_n_id + n_mc
		extend_nodes_t_front[[t]] = next_n_id:(next_n_id + n_mc - 1)
		next_n_id = next_n_id + n_mc
		growth_nodes_t[[t]] = next_n_id
		next_n_id = next_n_id + 1
	}
#we eliminate the last growth node by overiding its id with the sink
	sink_id = next_n_id - 1

#add edges from source to back time 1
	next_e_id = 1
	edges_src = data.frame(
						from = rep(source_id, n_mc),
						to = normal_nodes_t_back[[1]],
						ID = next_e_id:(next_e_id + n_mc - 1),
						capacity = rep(k_inf, n_mc),
						min_capacity = rep(0,n_mc),
						cost = rep(0, n_mc),
						mc1 = rep(-1, n_mc),
						mc2 = 1:n_mc,
						time1 = rep(0, n_mc), time2 = rep(1, n_mc),
						type1 = rep("src", n_mc), type2 = rep("norm_b", n_mc))
	next_e_id = next_e_id + n_mc

	edges_capacity = list()
	edges_manifold = list()
	edges_growth_in = list()
	edges_growth_out = list()
	alpha = capacity_var_factor
	for(t in 1:max_t) {
		edges_capacity[[t]] = data.frame(
						from = c(rep(normal_nodes_t_back[[t]],2), 
										rep(extend_nodes_t_back[[t]],2)),
						to = c(rep(normal_nodes_t_front[[t]],2), 
										rep(extend_nodes_t_front[[t]],2)),
						ID = next_e_id:(next_e_id + 4*n_mc - 1),
						capacity = c((1-1/(1+alpha))*mc_t_freq[,t],
										 (1/(1+alpha))*mc_t_freq[,t],
										 (alpha)*mc_t_freq[,t],
										 ifelse(mc_t_freq[,t]>0, k_inf,0)),
						min_capacity = rep(0, 4*n_mc),
						cost = rep(c(-off_capacity_cost2, -off_capacity_cost1, 
										off_capacity_cost1, off_capacity_cost2), each=n_mc),
						mc1 = rep(1:n_mc, 4),
						mc2 = rep(1:n_mc, 4),
						time1 = rep(t, 4*n_mc), time2 = rep(t, 4*n_mc),
						type1 = c(rep("norm_b", 2*n_mc), rep("extend_b", 2*n_mc)),
						type2 = c(rep("norm_f", 2*n_mc), rep("extend_f", 2*n_mc)))

		next_e_id = next_e_id + 4*n_mc

		if(t != max_t) {
			n_e = nrow(manifold)
			norm_norm_from = normal_nodes_t_front[[t]][manifold$mc1]
			norm_norm_to = normal_nodes_t_back[[t+1]][manifold$mc2]

			norm_ext_from = normal_nodes_t_front[[t]][manifold$mc1]
			norm_ext_to = extend_nodes_t_back[[t+1]][manifold$mc2]

			ext_norm_from = extend_nodes_t_front[[t]][manifold$mc1]
			ext_norm_to = normal_nodes_t_back[[t+1]][manifold$mc2]

			ext_ext_from = extend_nodes_t_front[[t]][manifold$mc1]
			ext_ext_to = extend_nodes_t_back[[t+1]][manifold$mc2]

			edges_manifold[[t]] = data.frame(
						from = c(norm_norm_from, norm_ext_from, ext_norm_from, ext_ext_from),
						to = c(norm_norm_to, norm_ext_to, ext_norm_to, ext_ext_to),
						ID = next_e_id:(next_e_id+4*n_e-1),
						capacity = rep(k_inf, 4*n_e),
						min_capacity = rep(0,4*n_e),
						cost = c(manifold$cost, 
									manifold$cost*k_norm_ext_cost, 
									manifold$cost*k_ext_norm_cost,
									manifold$cost*k_ext_ext_cost),
						mc1 = rep(manifold$mc1, 4),
						mc2 = rep(manifold$mc2, 4),
						time1 = rep(t, 4*n_e),
						time2 = rep(t+1, 4*n_e),
						type1 = c(rep("norm_f", 2*n_e), rep("extend_f", 2*n_e)),
						type2 = rep(c(rep("norm_b", n_e), rep("extend_b", n_e)),2))

			next_e_id = next_e_id + 4*n_e

			edges_growth_in[[t]] = data.frame(
						from = normal_nodes_t_front[[t]],
						to = rep(growth_nodes_t[[t]], n_mc),
						ID = next_e_id:(next_e_id + n_mc - 1),
						capacity = growth_leak[, t],
						min_capacity = rep(0,n_mc),
						cost = rep(0, n_mc),
						mc1 = 1:n_mc,
						mc2 = -2,
						time1 = rep(t, n_mc), time2 = rep(t, n_mc),
						type1 = rep("norm_f", n_mc), type2 = rep("growth", n_mc))
			next_e_id = next_e_id + n_mc
			edges_growth_out[[t]] = data.frame(
						from = rep(growth_nodes_t[[t]], n_mc),
						to = normal_nodes_t_back[[t+1]],
						ID = next_e_id:(next_e_id + n_mc - 1),
						capacity = growth_compensate[,t+1],
						min_capacity = rep(0,n_mc),
						cost = rep(0, n_mc),
						mc1 = -2,
						mc2 = 1:n_mc,
						time1 = rep(t, n_mc), time2 = rep(t+1, n_mc),
						type1 = rep("growth", n_mc), type2 = rep("norm_b", n_mc))
			next_e_id = next_e_id + n_mc
		}	
	}
	edges_sink = data.frame(
						from = c(normal_nodes_t_front[[max_t]],
									extend_nodes_t_front[[max_t]]),
						to = rep(sink_id, n_mc*2),
						ID = next_e_id:(next_e_id + n_mc*2 - 1),
						capacity = rep(k_inf, 2*n_mc),
						min_capacity = rep(0, 2*n_mc),
						cost = rep(0, 2*n_mc),
						mc1 = rep(1:n_mc, 2),
						mc2 = -1,
						time1 = rep(max_t, 2*n_mc), time2 = rep(max_t+1, 2*n_mc),
						type1 = c(rep("norm_f", n_mc), rep("extend_f", n_mc)), type2 = rep("sink", 2*n_mc))
	next_e_id = next_e_id + n_mc

	ec = do.call("rbind", edges_capacity)
	em = do.call("rbind", edges_manifold)
	egi = do.call("rbind", edges_growth_in)
	ego = do.call("rbind", edges_growth_out)

	network = rbind(edges_src, ec)
	network = rbind(network, em)
	network = rbind(network, egi)
	network = rbind(network, ego)
	network = rbind(network, edges_sink)
	
	
	# construct dataframe with all node IDs
	t_node = c()
	for (i in 1:max_t) {
	  t_node = c(t_node,rep(i,n_mc))
	}
	
	df_norm_ext_nodes = data.frame(node = c(unlist(normal_nodes_t_back),unlist(normal_nodes_t_front),unlist(extend_nodes_t_back),unlist(extend_nodes_t_front)),
	                               t = rep(t_node,4),
	                               mc = rep(c(1:n_mc),4*max_t),
	                               type = c(rep("norm_b",n_mc*max_t),rep("norm_f",n_mc*max_t),rep("extend_b",n_mc*max_t),rep("extend_f",n_mc*max_t)))
	
	growth_nodes = unlist(growth_nodes_t)[1:(max_t-1)]
	
	df_growth_nodes = data.frame(node = unlist(growth_nodes),
	                             t = c(1:length(growth_nodes)),
	                             mc = rep(-2,length(growth_nodes)),
	                             type = rep("growth",length(growth_nodes)))
	
	df_src_sink_nodes = data.frame(node = c(source_id,sink_id),t = c(0,max_t +1),mc = c(-1,-1),type = c("src","sink"))
	
	nodes = rbind(df_src_sink_nodes,df_growth_nodes)
	nodes = rbind(nodes,df_norm_ext_nodes)
	
	nodes = nodes[order(nodes$node),]
   rownames(nodes) = c(1:nrow(nodes))

	mct@network = network
	
#	return(list(net = network,
#	            nodes = nodes,
#					    norm_back_t_mc = normal_nodes_t_back,
#					    norm_front_t_mc = normal_nodes_t_front,
#					    ext_back_t_mc = extend_nodes_t_back,
#					    ext_front_t_mc = extend_nodes_t_front))

	return(mct)
}

build_growth_mats = function(mct, mc_leak)
{
	max_mc = nrow(mct@mc_t)
	min_t = 1
	max_t = ncol(mct@mc_t)
	mc_t_freq = mct@mc_t

	growth_leak = matrix(0, nrow=max_mc, ncol=(max_t-min_t+1))
	for(t in min_t:max_t) {
		growth_leak[,t-(min_t-1)] = mc_t_freq[,t]*mc_leak
	}

	tot_mass = colSums(mc_t_freq)

	lost_mass = colSums(growth_leak)/tot_mass

	growth_compensate = t(t(mc_t_freq)*c(0,lost_mass[-max_t]))
	return(list(leak=growth_leak, compensate=growth_compensate))
}

#' Compute mincost flow and update the mct objects with flows and 
#' inferred mc capacity per time 
#'
#' @param mct mct network object
#' @param flow_tolerance how much flow we should miss per time
#'
#' @return an updated mct object
#' @export
mctnetwork_gen_mincost_flows= function(mct, flow_tolerance=0.01)
{
	net = mct@network
	ncnstr = network_lp_constraints(net, 1-flow_tolerance, 100)

	edge_costs = net$cost[order(net$ID)]
	sol = lpsymphony::lpsymphony_solve_LP(obj = edge_costs,
										mat = ncnstr$lhs, 
										rhs = ncnstr$rhs, 
										dir = ncnstr$dir)
  
	flows = sol$solution
	net$flow = flows[net$ID]
	mct@network = net
	mct@edge_flows = net$flow
	names(mct@edge_flows) = net$ID	

	f_cap = net$mc1 == net$mc2 & net$time1 == net$time2
	net_cap = net[f_cap,]

	mc_time_weight = summarize(group_by(net_cap,mc1,time1),tot_weight = sum(flow))
	mc_t_post = sparseMatrix(i = as.numeric(mc_time_weight$mc1),
									j = as.numeric(mc_time_weight$time1), 
									x = mc_time_weight$tot_weight)
	mc_t_post = as.matrix(mc_t_post)
	rownames(mc_t_post) = c(1:nrow(mc_t_post))
	colnames(mc_t_post) = c(1:ncol(mc_t_post))

	
	mct@mc_t_infer = mc_t_post
	return(mct)
}


#' Compute matrices for forward and backwrd propagation using the flows
#'
#' @param mct mct network object
#'
#' @return an updated mct object
#' @export

mctnetwork_comp_propagation= function(mct)
{
#flows mc,mc,t
	max_t = ncol(mct@mc_t)
	for(t in 1:(max_t-1)) {
		mc_flows = mctnetwork_get_flow_mat(mct, t)

		mc_forward = mc_flows/rowSums(mc_flows)
		mc_backward= t(t(mc_flows)/colSums(mc_flows))
		mc_backward[is.nan(mc_backward)] = 0
		mc_forward[is.nan(mc_forward)] = 0
		mct@mc_forward[[t]] = mc_forward
		mct@mc_backward[[t]] = mc_backward
	}
	return(mct)
#forward flows
}

   
#' Compute matrix of flows between MCs in a given time point
#'
#' @param mct mct network object
#' @param time flows will be computed for the (time,time+1) interval. Time=-1 will generate total flow over all times
#'
#' @return a matrix on metacells
#' @export

mctnetwork_get_flow_mat = function(mct, time, max_time=time)
{
	net = mct@network
	if(time == -2) {
		f_t = net$type1 != "growth" & net$type2!="growth" & net$type2 != "sink"
	} else {
		if(time == -1) {
			f_t = net$type1 != "growth" & net$type2!="growth" & net$type1 != "source" & net$type2 != "sink"
		} else {
			f_t = net$time1 == time & net$time2==time+1 &
					net$type1 != "growth" & net$type2!="growth"
		}
	}

	net_t = net[f_t,] 
   flow = as.data.frame(summarize(group_by(net_t, mc1, mc2),
													tot_flow = sum(flow)))
    
   mc_mat = pivot_wider(data = flow, 
				names_from = mc2, 
				values_from = tot_flow,
				values_fill = list(tot_flow = 0))

   mc_mat = as.data.frame(mc_mat)
   rownames(mc_mat) = mc_mat$mc1
   mc_mat = mc_mat[,-1]
	max_mc = ncol(mc_mat)
	if(time == -2) {
		mc_mat = mc_mat[as.character(c(-1, 1:max_mc)), as.character(1:max_mc)]
		mc_mat = cbind(rep(0,nrow(mc_mat)), mc_mat)
	} else {
		mc_mat = mc_mat[as.character(1:max_mc), as.character(1:max_mc)]
	}
   mc_mat = as.matrix(mc_mat)

	return(mc_mat)
}

#' Compute backward and forward flow propgation of metacell probability from time t
#'
#' @param mct mct network object
#' @param t flows will be computed for the (time,time+1) interval
#' @param mc_p probabilities at time t
#'
#' @return a list with two elements: probs is a matrix of probabilities over metacells (rows)  and time (columns). step_m is a list of sparse matrices inferred flows between metacells per time.
#' @export

mctnetwork_propogate_from_t = function(mct, t, mc_p)
{
	max_t = ncol(mct@mc_t)
	step_m = list()
	probs = matrix(0, nrow = nrow(mct@mc_t), ncol=ncol(mct@mc_t))
	probs[,t] = mc_p
	if(t > 1) {
		for(i in (t-1):1) {
			step_m[[i]] = Matrix(t(t(as.matrix(mct@mc_backward[[i]])) * probs[,i+1]), sparse=T)
			
			probs[,i] = as.matrix(mct@mc_backward[[i]]) %*% probs[,i+1]
		}
	}
	if(t < max_t) {
		for(i in (t+1):max_t) {
			step_m[[i-1]] = Matrix(as.matrix(mct@mc_forward[[i-1]]) * probs[,i-1], sparse=T)
			probs[,i] = t(probs[,i-1]) %*% as.matrix(mct@mc_forward[[i-1]])
		}
	}
	return(list(probs=probs, step_m=step_m))
}

#' Given a list of genes and a matrix, this compute mean umi per metacell (rows) and time (column).
#'
#' @param mct_id  mct network object
#' @param mat_id  umi matrix object
#' @param genes list of gene names 
#' @param min_percentile percentile value to use as minimum threshold (this is used to avoid very small values and regularize color scale, default 0.05)
#'
#' @return A matrix on metacells and time with mean umi fraction over the gene module
#' @export
#'
mctnetwork_gen_gmod_mc_t = function(mct_id, mat_id,  genes, min_percentile=0.05)
{
	mct = scdb_mctnetwork(mct_id)
	if(is.null(mct)) {
		stop("cannot find mctnet object ", mct_id, " when trying to plot net flows anchors")
	}
	mc = scdb_mc(mct@mc_id)
	mat = scdb_mat(mat_id)
	genes = intersect(genes, rownames(mat@mat))
	gmod_tot = colSums(mat@mat[genes,names(mc@mc)])/colSums(mat@mat[,names(mc@mc)])
	gmod_df = data.frame(mc = mc@mc, 
								t = mct@cell_time[names(mc@mc)], 
								gmod = gmod_tot)
	mc_time_gmod = summarise(group_by(gmod_df, mc, t), gmod=mean(gmod))
	mc_t_gmod_m = sparseMatrix(i=mc_time_gmod$mc, 
									j=mc_time_gmod$t, 
									x=mc_time_gmod$gmod)
	mc_t_gmod_m = pmax(mc_t_gmod_m, quantile(mc_time_gmod$gmod, min_percentile))
	return(mc_t_gmod_m)
}

#' This generate two large matrices showing expression of genes per metacell ad time point, as well as the expression of the gene in inferred ancestral states given the flow model
#'
#' @param mct_id  mct network object
#' @param mat_id  umi matrix object
#' @param genes list of gene names 
#'
#' @return A matrix on metacells and time with mean umi fraction over the gene module
#' @export
#'
mctnetwork_gen_e_gmt_p = function(mct_id, mat_id,  genes)
{
	mct = scdb_mctnetwork(mct_id)
	if(is.null(mct)) {
		stop("cannot find mctnet object ", mct_id, " when trying to plot net flows anchors")
	}
	mc = scdb_mc(mct@mc_id)
	mat = scdb_mat(mat_id)
	genes = intersect(genes, rownames(mat@mat))

	max_t = length(mct@mc_backward) + 1

	csize = colSums(mat@mat)
	e_gm_t = list()
	tot_m_t = list()
	for(t in 1:max_t) {
		cell_t = intersect(names(mc@mc), names(which(mct@cell_time==t)))
		message("at ", t, " with ", length(cell_t), " cells")
		mgt = as.matrix(mat@mat[genes, cell_t])
#		tot_gm = t(tgs_matrix_tapply(mgt, mc@mc[cell_t], sum))
		tot_gm = t(apply(mgt,1, function(x) tapply(x,mc@mc[cell_t], sum)))
		tot_c = tapply(csize[cell_t], mc@mc[cell_t], sum)
		tot_m_t[[t]] = tot_c
		e_gm_t[[t]] = t(t(tot_gm)/as.vector(tot_c))
		rownames(e_gm_t[[t]]) = genes
	}
	e_gm_t_prev = list()
	for(t in 2:max_t) {
		t_mcs = colnames(e_gm_t[[t-1]])
		t_mcs2= colnames(e_gm_t[[t]])
		e_gm_t_prev[[t]] = Matrix(e_gm_t[[t-1]] %*% as.matrix(mct@mc_backward[[t-1]][t_mcs, t_mcs2]), sparse=T)
	}
	
	e_gm_t_prev2 = list()
	for(t in 3:max_t) {
	  t_mcs0 = colnames(e_gm_t[[t-2]])
	  t_mcs = colnames(e_gm_t[[t-1]])
	  t_mcs2= colnames(e_gm_t[[t]])
	  e_gm_t_prev2[[t]] = Matrix(e_gm_t[[t-2]] %*% as.matrix(mct@mc_backward[[t-2]][t_mcs0, t_mcs]) %*% as.matrix(mct@mc_backward[[t-1]][t_mcs, t_mcs2]), sparse=T)
	}
	
	return(list(e_gm_t = e_gm_t, e_gm_t_prev = e_gm_t_prev, e_gm_t_prev2 = e_gm_t_prev2, tot_m_t = tot_m_t))
}

#' Compute matrix of flows over cell types
#'
#' @param mct mct network object
#' @param min_time minimum time point
#' @param max_time maximum time point
#'
#' @return a list of matrices show flows from type t to t'
#' @export

mctnetwork_get_type_flows = function(mct, time, max_time)
{
	net = mct@network
	mc = scdb_mc(mct@mc_id)

	all_types = unique(mc@colors)
	mct_mats = list()
	for(t in time:(max_time-1)) {
		f_t = net$time1 == t & net$time2 == t+1 &
					net$type1 != "growth" & net$type2!="growth"

		net_t = net[f_t,] 
		net_t$mc_t1 = mc@colors[as.numeric(net_t$mc1)]
		net_t$mc_t2 = mc@colors[as.numeric(net_t$mc2)]
		flow = as.data.frame(summarize(group_by(net_t, mc_t1, mc_t2),
													tot_flow = sum(flow)))
    
	   mct_mat = pivot_wider(data = flow, 
				names_from = mc_t2, 
				values_from = tot_flow,
				values_fill = list(tot_flow = 0))

	   mct_mat = as.data.frame(mct_mat)
		rownames(mct_mat) = mct_mat$mc_t1
	   mct_mat = mct_mat[,-1]
		max_mc = ncol(mct_mat)
		mct_mat = mct_mat[all_types, all_types]
		mct_mat = as.matrix(mct_mat)
		mct_mats[[t]] = mct_mat
	}

	return(mct_mats)
}

mctnetwork_get_egc_on_cluster_transition = function(mct, min_time, max_time, type1, type2, mc_type=NULL)
{
	mc = scdb_mc(mct@mc_id)
	e_gc = mc@e_gc
	net = mct@network

	if(is.null(mc_type)) {
		mc_type = mc@colors
		names(mc_type) = as.character(1:length(mc_type))
	}

#	flow_mm = mctnetwork_get_flow_mat(mct, time, max_time=time)

	f_t = net$time1 >= min_time & net$time2 <= max_time &
					net$time1 == net$time2-1 &
					net$type1 != "growth" & net$type2!="growth"

	net = net[f_t,]
	f_types = mc_type[as.numeric(net$mc1)]==type1 & mc_type[as.numeric(net$mc2)]==type2

	net = net[f_types,]

	src_mc_wgt = tapply(net$flow, net$mc1, sum)	
	targ_mc_wgt = tapply(net$flow, net$mc2, sum)
	src_mc_wgt_n = as.vector(src_mc_wgt/sum(src_mc_wgt))
	names(src_mc_wgt_n) = names(src_mc_wgt)
	targ_mc_wgt_n = as.vector(targ_mc_wgt/sum(targ_mc_wgt))
	names(targ_mc_wgt_n) = names(targ_mc_wgt)

	src_e_gc = colSums(t(e_gc[,names(src_mc_wgt_n)]) * src_mc_wgt_n)
	targ_e_gc = colSums(t(e_gc[,names(targ_mc_wgt_n)]) * targ_mc_wgt_n)

	return(data.frame(src = src_e_gc, targ = targ_e_gc, lf = log2(1e-5+targ_e_gc)-log2(1e-5+src_e_gc)))
}
