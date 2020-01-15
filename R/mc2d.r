#' Meta cell 2d projection
#'
#' Representing a graph of metacells and their 2D embeddings
#'
#' @slot mc_id id of the metacell object
#' @slot x x coordinate of a (subset) of the metacells
#' @slot y y coordinate of a (subset) of the metacells
#' @slot graph datafame with the edges linking metacells
#'
#'
#' @export tgMC2D
#' @exportClass tgMC2D

tgMC2D <- setClass(
   "tgMC2D",
	slots = c(
	  mc_id = "character",
	  mc_x = "vector",
	  mc_y = "vector",
	  sc_x = "vector",
	  sc_y = "vector",
	  graph = "data.frame")
)

#' Construct a meta cell 2D embedding
#'
#' Just filling up the slots, actuall computation is done by specific algorithms (mcell_mc2d_force)
#'
#' @param mc metacell object
#' @param mc_x x coordinates of metacells
#' @param mc_y y coordinates of metacells
#' @param sc_x x coordinates of cells
#' @param sc_y y coordinates of cells
#' @param graph data.frame with fields mc1, mc2 definining the graph
#' @export

setMethod(
  "initialize",
  signature = "tgMC2D",
  definition =
    function(.Object, mc_id, mc_x, mc_y, sc_x, sc_y, graph) {
		if(is.null(scdb_mc(mc_id))) {
			stop("initlaizing mc2d using mc id ", mc_id, " that is missing from scdb")
		}
		.Object@mc_id = mc_id
		.Object@mc_x = mc_x
		.Object@mc_y = mc_y
		.Object@sc_x = sc_x
		.Object@sc_y = sc_y
		if(length(intersect(colnames(graph), c("mc1", "mc2"))) != 2) {
			stop("Bad tgMC2D initialization, graph parameter must be a dta frame with fields mc1, mc2")
		}
		.Object@graph= graph
      return(.Object)
    }
)

#' Generate a new metacell in scdb
#'
#' This constructs a meta cell cover object and puts it into scdb. It gets an MC assignment (cell->MC_ID), and a matrix, and call standard api of this class to compute the footprints.
#'
#' @param mc2d_id new id of 2d object to add
#' @param mc_id id of scdb meta cell object we link to
#' @param mc_x x coordinates of metacells
#' @param mc_y y coordinates of metacells
#' @param sc_x x coordinates of cells
#' @param sc_y y coordinates of cells
#' @param graph graph data frame (fields mc1,mc2 specifying edges between mc's)
#' @export
mcell_new_mc2d = function(mc2d_id, mc_id, mc_x, mc_y, sc_x, sc_y, graph)
{
	scdb_add_mc2d(mc2d_id, tgMC2D(mc_id, mc_x, mc_y, sc_x, sc_y, graph))
}

#' Rotatae/invert projection
#'
#' transform all coordainte by a rotation and inversion
#'
#' @param mc2d_id new id of 2d object to add
#' @param alpha angle in degrees
#' @param flipx should flip the X over the mean(x)
#' @param flipy should flip the Y over the mean(Y)
#'
#' @export
mcell_mc2d_rotate= function(mc2d_id, alpha = 0, flipx = F, flipy = F)
{
	mc2d = scdb_mc2d(mc2d_id)
	if(is.null(mc2d)) {
		stop("cannot get mc2d object id ", mc2d_id)
	}
	alpha = alpha*pi/180
	sina = sin(alpha)
	cosa = cos(alpha)
	mx = mean(mc2d@sc_x, na.rm=T)
	my = mean(mc2d@sc_y, na.rm=T)
	if(flipx) {
		mc2d@sc_x = 2*mx - mc2d@sc_x
		mc2d@mc_x = 2*mx - mc2d@mc_x
	}
	if(flipy) {
		mc2d@sc_y = 2*my - mc2d@sc_y
		mc2d@mc_y = 2*my - mc2d@mc_y
	}

	dscx = mc2d@sc_x - mx
	dscy = mc2d@sc_y - my
	scx = mx + cosa*dscx - sina*dscy
	scy = my + sina*dscx + cosa*dscy

	dmcx = mc2d@mc_x - mx
	dmcy = mc2d@mc_y - my
	mcx = mx + cosa*dmcx - sina*dmcy
	mcy = my + sina*dmcx + cosa*dmcy

	mc2d@sc_x = scx
	mc2d@sc_y = scy
	mc2d@mc_x = mcx
	mc2d@mc_y = mcy

	scdb_add_mc2d(mc2d_id, mc2d)
}

