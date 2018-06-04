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

#' Plot a color bar with values
#'
#' @param vals to plot
#' @param cols colors of the values (same length as vals)
#' @param fig_fn if null - plot inline
#' @param title of legend
#' @param show_vals_ind logical vector - which values to show
#'
plot_color_bar = function(vals, cols, fig_fn=NULL, title="", show_vals_ind=NULL)
{
  if (!is.null(fig_fn)) {
    .plot_start(fig_fn, 400, 400)
  }
  plot.new()
  plot.window(xlim=c(0,100), ylim=c(0, length(cols) + 3))
  rect(7, 1:length(cols), 17, 1:length(cols) + 1, border=NA, col=cols)
  rect(7, 1, 17, length(cols)+1, col=NA, border = 'black')

  if (is.null(show_vals_ind)) {
    show_vals_ind = rep(T, length(cols))
  }
  text(19, (1:length(cols))[show_vals_ind] + 0.5, labels=vals[show_vals_ind], pos=4)
  text(2, length(cols)/2 + 1, labels=title, srt=90, cex=1.5)

  if (!is.null(fig_fn)) {
    dev.off()
  }
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
