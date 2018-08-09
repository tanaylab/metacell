.onLoad <- function(libname, pkgname) {
	tgconfig::register_params(system.file('config/metacell_params.yaml', package='metacell'), package='metacell')
	ggplot2::theme_set(ggplot2::theme_classic())
}
