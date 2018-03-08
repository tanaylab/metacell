#' scgfig_set_base - set base directory for metacell figures, creates it if it doesn't exist.
#'
#' @export

scfigs_init = function(base)
{
	dir.create(base, recursive = T, showWarnings = F)
	.scfigs_base <<- base
}

#' Generate a standard figure name igven and object and figure type
#'
#' @param id - is of the object the figure relates to (E.g. mat_id)
#' @param type - a string defining the figure type
#' @param dir - output dir. If null, using figs base dir. Creates it if it doesn't exist.
#'
#' @export
#'
scfigs_fn = function(id, type, dir = NULL)
{
  if (is.null(dir)) {
    dir = .scfigs_base
  }

	if(!file.exists(dir)) {
		dir.create(dir, recursive = T, showWarnings = F)
	}

	return(sprintf("%s/%s.%s.png", dir, id, type))
}

#' Generate a standard figure dir name igven and object and figure type
#'
#' @param id - is of the object the figure relates to (E.g. mat_id)
#' @param type - a string defining the figure type
#'
#' @export
#'
scfigs_dir = function(id, type)
{
	if(is.null(.scfigs_base) | !file.exists(.scfigs_base)) {
		stop("figs directory at ", .scfigs_base, " is missing")
	}
	dir_nm = sprintf("%s/%s.%s", .scfigs_base, id, type)
	if(!dir.exists(dir_nm)) {
		dir.create(dir_nm)
	}
	return(dir_nm)
}

# wrap for opening a plot (png, ps or pdf)
.plot_start = function(fn, w, h)

{
  device = get_param("mc_plot_device")
  res = get_param("mc_plot_ppi")

	if (device == "png") {
		png(filename=sub("ps$", "png", fn), width=w, height=h, res=res)
	}
	else if (device == "ps") {
		postscript(file=sub("png$", "ps", fn), width=w/res, height=h/res)
	}
	else if (device == "pdf") {
		pdf(file=sub("png$", "pdf", fn), width=w/res, height=h/res)
	}
	else {
		stop(sprintf("unknown output device type: %s", device))
	}
}
