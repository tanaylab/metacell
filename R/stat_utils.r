

#' port_cor
#' matrix correlation wrap. Parameter selects which function to use, cor or tgs_cor. Tgs_cor can be more efficien if available
#'
#' @param x matrix whole columans will be analyzed
#' @param pairwise.complete.obs - similar to cor, if this is T, NA's will be ibe ignored and only values that are pairwise available will be used to compute the correlaiton
#' @param spearman TRUE if rank correlation is to be used
#' @param tidy if TRUE, the output will be tidy
#' @param threshold defining a threshold to deterine which correlation to be reported if using the tidy format
#'
#' @return the corelation matrix of the columns
#'
port_cor = function(x, pairwise.complete.obs = F, spearman = F, tidy = F,	threshold = 0)
{
	if (get_param("mc_use_tgs_cor")) {
		return (tgs_cor(x, pairwise.complete.obs, spearman, tidy, threshold))
	}
	else {
		if(tidy) {
			stop("MC-ERR no tidy support for correaltion when not using tgs_cor, try changing the param mc_use_tgs_cor to TRUE")
		}
		return(cor(x, method=ifelse(spearman, "spearman", "pearson"), use=ifelse(pairwise.complete.obs, "pair", "all")))
	}
}
