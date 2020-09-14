
#edges should be a data frame with
#from
#to,
#ID
#capacity
#cost

network_lp_constraints <- function(net, total_flow, k_inf_cap) 
{

  f_capped_edges = net$capacity < k_inf_cap
  n_cap = sum(f_capped_edges)
  lhs_cap_constr = data.frame(id=1:n_cap,
									edge= net$ID[f_capped_edges],
									coef=rep(1, n_cap))
  
  rhs_cap_constr = net$capacity[f_capped_edges]
  dir_cap_constr = rep('<=', times = n_cap)

  next_cnstr_id = nrow(lhs_cap_constr) + 1
  
  f_min_capped_edges = net$min_capacity > 0
  n_min_cap = sum(f_min_capped_edges)
  if(n_min_cap > 0) {
		lhs_min_cap_constr = data.frame(id = next_cnstr_id:(next_cnstr_id + n_min_cap - 1),
                                  edge = net$ID[f_min_capped_edges],
                                  coef = rep(1,n_min_cap))
  
		rhs_min_cap_constr = net$min_capacity[f_min_capped_edges]
		dir_min_cap_constr = rep('>=', times = n_min_cap)
  
		next_cnstr_id = next_cnstr_id + nrow(lhs_min_cap_constr)
		lhs_cap_constr = rbind(lhs_cap_constr, lhs_min_cap_constr)
		rhs_cap_constr = c(rhs_cap_constr, rhs_min_cap_constr)
		dir_cap_constr = c(dir_cap_constr, dir_min_cap_constr)
  }
  
#no source constraint
  f = net$from != 1
  lhs_flow_constr_o = data.frame(id = (next_cnstr_id - 2 + net$from[f]),
									edge = net$ID[f],
									coef = rep(-1,sum(f))) 

#no sink consraint
  f = net$to != max(net$to)
  lhs_flow_constr_i = data.frame(id =next_cnstr_id -2 + net$to[f],
									edge = net$ID[f],
									coef=rep(1, sum(f)))
  rhs_flow_constr = rep(0, max(net$to)- 2)
  dir_flow_constr = rep('==', max(net$to) - 2)

  next_cnstr_id = next_cnstr_id + max(net$to) - 2

  src_adj = net[net$from == 1,"ID"]
  sink_adj = net[net$to == max(net$to),"ID"]
  lhs_src_constr = data.frame(id = rep(next_cnstr_id, length(src_adj)),
									edge = src_adj,
									coef=rep(1, length(src_adj)))
  rhs_src_constr = total_flow
  dir_src_constr = '=='

  next_cnstr_id = next_cnstr_id + 1
  lhs_sink_constr = data.frame(id = rep(next_cnstr_id, length(sink_adj)),
									edge = sink_adj,
									coef=rep(1, length(sink_adj)))
  rhs_sink_constr = total_flow
  dir_sink_constr = '=='

  # Build constraints matrix
  lhs_df = rbind(rbind(lhs_cap_constr, lhs_flow_constr_o), lhs_flow_constr_i)
  lhs_df = rbind(rbind(lhs_df, lhs_src_constr), lhs_sink_constr)
#  lhs = sparseMatrix(i = lhs_df$id, lhs_df$edge, x= lhs_df$coef)
  lhs = slam::simple_triplet_matrix(i = lhs_df$id, lhs_df$edge, v= lhs_df$coef)
  constraints <- list(
    lhs = lhs,
    dir = c(dir_cap_constr, dir_flow_constr, dir_src_constr, dir_sink_constr),
    rhs = c(rhs_cap_constr, rhs_flow_constr, rhs_src_constr, rhs_sink_constr))
	return(constraints)
}

gen_lp_from_net = function(net, ncnstr) 
{
#	lpfl = make.lp(nrow=nrow(ncnstr$lhs), ncol=ncol(ncnstr$lhs))
	lpfl = make.lp(nrow=nrow(ncnstr$lhs))
	for(i in seq(1,ncol(ncnstr$lhs),10000)) {
		to_i = min(i + 10000 - 1, ncol(ncnstr$lhs))
		message("add ", i, " - ", to_i, " cnstr")
		apply(ncnstr$lhs[,i:to_i], 2, function(x) add.column(lpfl, x))
	}
	set.objfn(lpfl, net$net$cost)
	set.constr.type(lpfl, ncnstr$dir)
	set.rhs(lpfl, ncnstr$rhs)

#	for(i in 1:nrow(ncnstr$lhs)) { 
#			add.constraint(lpfl, ncnstr$lhs[i,], ncnstr$dir[i], ncnstr$rhs[i]) 
#	}
	opt = lp.control(lpfl, verbose="normal")
	return(lpfl)
}
set_lp_net_maxflow = function(lpfl, maxflow)
{
	sink_id = length(get.constr.type(lpfl))
	src_id = length(get.constr.type(lpfl)) - 1
	set.constr.value(sink_id, rsh = maxflow, constraints = c(sink_id, src_id))
}

