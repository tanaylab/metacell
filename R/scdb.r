#' Initializing scdb
#'
#' This initialize scdb to a certain directory. The object will then allow loading. chaching and saving of different objects, including matrix, gene statictics, cell graphs, metacell covers, and more
#' 
#'
#' @examples
#' # we first initialize a db object
#' scdb_init("workdir")
#  # we can see which objects are available
#' scdb_ls("mat")
#' TBA
#'
#' @export

scdb_init = function(base_dir, force_reinit=F)
{
	if(exists(".scdb") & !force_reinit) {
		stop("scdb already init to ", 
				ifelse(exists(".scdb_base"), .scdb_base, "NA"), 
				" use force reinit to restart it")
	}
	if(!dir.exists(base_dir)) {
		stop("cannot initialize db to non existing directory base_dir")
	} else {
		message("initializing scdb to ", base_dir)
	}
	.scdb_base <<- base_dir
#init repos
	.scdb <<- list("mat" = list(), 
		     "gstat" = list(),
		     "gset" = list(),
		     "mc" = list(), 
		     "cgraph" = list(), 
		     "coclust" = list(),
		     "mc2d" = list())
}

#' Testing is scdb is initialized
#'
#' Retruning TRUE if scdb is initialized
#'
#' @export
#'
scdb_is_valid = function()
{
	if(!exists(".scdb")) {
		return(FALSE)
	} else {
		if(file.exists(.scdb_base)) {
			return(TRUE)
		} else {
			stop(".scdb ", .scdb, " points to non existing directory")
		}
	}
}

#' List all object of a given type from the current scdb
#'
#' @param objt - the object type
#' @export
#'
scdb_ls = function(objt)
{
	if(!exists(".scdb")) {
		message("scdb not initialized")
	} else {
		fns = grep(sprintf("^%s",objt),
					list.files(sprintf("%s/", .scdb_base)), 
					v=T)
		fns = sub(".Rda","", fns)
		print(fns)
	}
}

#' scdb_ls_loaded - list loaded object of a certain type
#'
#' @param objt - either mat, 
#'
#'
#' @export

scdb_ls_loaded = function(objt)
{
	if(!exists(".scdb")) {
		message("scdb not initialized")
	} else {
		print(names(.scdb[[objt]]))
	}
}

.scdb_obj_fn = function(objt, id)
{
	fn = sprintf("%s/%s.%s.Rda", .scdb_base, id, objt)
	return(fn)
}

.scdb_get_obj = function(id, objt) 
{
	if(!exists(".scdb")) {
		message("scdb not initialized, cannot get object ", id, " type ", objt)
		return(NULL)
	}
	if(is.null(.scdb[[objt]][[id]])) {
		fn = .scdb_obj_fn(id, objt)
		if(file.exists(fn)) {
			load(fn)
			.scdb[[objt]][[id]] <<- object
			return(object)
		} else {
			return(NULL)
		}
	} else {
		return(.scdb[[objt]][[id]])
	}
}

.scdb_add_obj = function(id, objt, object)
{
	.scdb[[objt]][[id]] <<- object
	fn = .scdb_obj_fn(id, objt)
	save(object, file=fn)
}

.scdb_del_obj = function(id, objt)
{
	.scdb[[objt]][[id]] <<- NULL
	fn = .scdb_obj_fn(id, objt)
	file.remove(fn)
}

#' scdb_mat - get matrix from db, load it if needed
#'
#' @param id -  matrix id. Will return null if does not exist
#'
#' @export

scdb_mat= function(id) 
{
	return(.scdb_get_obj(id, "mat"));
}

#' scdb_add_mat - add amatrix to the DB - will save it and cache
#'
#' @param id -  matrix id. Will return null if does not exist
#' @param mat  - matrix object
#'
#' @export
#'
scdb_add_mat= function(id, mat) 
{
#	if(typeof(mat) != "tgScMat") {
#		stop("Cannot add non tgScMat object as a mat in scdb")
#	}
	.scdb_add_obj(id, "mat", mat);
}

#' scdb_del_mat - remove a matrix from the DB (not just the cache!)
#'
#' @param id -  matrix id to remove from the DB
#'
#' @export
#'
scdb_del_mat= function(id)
{
	.scdb_del_obj(id, "mat");
}

#' scdb_gstat - get a gstat data frame. If it is missing and the id match an existing matrix, a gstat will be gerated for this matrix and added to scdb
#'
#' @param id - id of gstat
#'
#' @export
#'
scdb_gstat = function(id) 
{
	return(.scdb_get_obj(id, "gstat"))
}

#' scdb_add_gstat - add gstat to the DB and cahce
#'
#' @param id - id of gstat
#' @param gstat - gstat data frame
#'
#' @export
#'
scdb_add_gstat = function(id, gstat) 
{
	if(typeof(gstat ) != "list") {
		stop("Cannot add non dataframe/list object as a gstat in scdb")
	}
	.scdb_add_obj(id, "gstat", gstat);
}

#' scdb_del_gstat - del gstat from the DB and cahce
#'
#' @param id - id of gstat
#'
#' @export
scdb_del_gstat = function(id)
{
	.scdb_del_obj(id, "gstat");
}

#' scdb_gset - get a gene set
#'
#' @param id - id of gset
#'
#' @export
#'
scdb_gset = function(id) 
{
	return(.scdb_get_obj(id, "gset"))
}

#' scdb_add_gset - add gset to the DB and cahce
#'
#' @param id - id of gset
#' @param gset - gset data frame
#'
#' @export
#'
scdb_add_gset = function(id, gset) 
{
	if(class(gset)[1] != "tgGeneSets") {
		stop("Cannot add non tgGeneSets object as a gset in scdb")
	}
	.scdb_add_obj(id, "gset", gset);
}

#' scdb_del_gset - del gset from the DB and cahce
#'
#' @param id - id of gset
#'
#' @export
scdb_del_gset = function(id)
{
	.scdb_del_obj(id, "gset");
}

#' scdb_cgraph - get a cgraph object
#'
#' @param id - id of cgraph
#'
#' @export
#'
scdb_cgraph = function(id) 
{
	return(.scdb_get_obj(id, "cgraph"))
}

#' scdb_add_cgraph - add cgraph to the DB and cahce
#'
#' @param id - id of cgraph
#' @param cgraph - cgraph data frame
#'
#' @export
#'
scdb_add_cgraph = function(id, cgraph) 
{
	if(class(cgraph)[1] != "tgCellGraph") {
		stop("Cannot add non tgCellGraph object as a cgraph in scdb")
	}
	.scdb_add_obj(id, "cgraph", cgraph);
}

#' scdb_del_cgraph - del cgraph from the DB and cahce
#'
#' @param id - id of cgraph
#'
#' @export
scdb_del_cgraph = function(id)
{
	.scdb_del_obj(id, "cgraph");
}

#' scdb_mc - get a mc object
#'
#' @param id - id of mc
#'
#' @export
#'
scdb_mc = function(id) 
{
	return(.scdb_get_obj(id, "mc"))
}

#' scdb_add_mc - add mc to the DB and cahce
#'
#' @param id - id of mc
#' @param mc - mc data frame
#'
#' @export
#'
scdb_add_mc = function(id, mc) 
{
	if(class(mc)[1] != "tgMCCov") {
		stop("Cannot add non tgMCCov object as a mc in scdb")
	}
	.scdb_add_obj(id, "mc", mc);
}

#' scdb_del_mc - del mc from the DB and cahce
#'
#' @param id - id of mc
#'
#' @export
scdb_del_mc = function(id)
{
	.scdb_del_obj(id, "mc");
}

#' scdb_mc2d - get a mc2d object
#'
#' @param id - id of mc2d
#'
#' @export
#'
scdb_mc2d = function(id) 
{
	return(.scdb_get_obj(id, "mc2d"))
}

#' scdb_add_mc2d - add mc2d to the DB and cahce
#'
#' @param id - id of mc2d
#' @param mc2d - mc2d data frame
#'
#' @export
#'
scdb_add_mc2d = function(id, mc2d) 
{
	if(class(mc2d)[1] != "tgMC2D") {
		stop("Cannot add non tgMC2D object as a mc2d in scdb")
	}
	.scdb_add_obj(id, "mc2d", mc2d);
}

#' scdb_del_mc2d - del mc2d from the DB and cahce
#'
#' @param id - id of mc2d
#'
#' @export
scdb_del_mc2d = function(id)
{
	.scdb_del_obj(id, "mc2d");
}

#' scdb_coclust - get a coclust object
#'
#' @param id - id of coclust
#'
#' @export
#'
scdb_coclust = function(id) 
{
	return(.scdb_get_obj(id, "coclust"))
}

#' scdb_add_coclust - add coclust to the DB and cahce
#'
#' @param id - id of coclust
#' @param coclust - coclust data frame
#'
#' @export
#'
scdb_add_coclust = function(id, coclust) 
{
	if(class(coclust)[1] != "tgCoClust") {
		stop("Cannot add non tgCoClust object as a coclust in scdb")
	}
	.scdb_add_obj(id, "coclust", coclust);
}

#' scdb_del_coclust - del coclust from the DB and cahce
#'
#' @param id - id of coclust
#'
#' @export
scdb_del_coclust = function(id)
{
	.scdb_del_obj(id, "coclust");
}
