#' MAE plot
#' @name gemini_plot_mae
#'
#' @description Plots the average MAE of gemini model at each iteration step.
#'
#' @param Model a gemini.model object
#'
#' @return a ggplot2 object
#' @import ggplot2
#' @importFrom scales pretty_breaks
#' @export
#' 
#' @examples
#' data("Model", package = "gemini")
#' gemini_plot_mae(Model)
#' 
gemini_plot_mae <- function(Model){
	stopifnot("gemini.model" %in% class(Model))
	stopifnot("mae" %in% names(Model))

	g = ggplot(data = data.frame(iter = seq_len(length(Model$mae))-1, mae = Model$mae),
					aes_string(x = 'iter', y = 'mae')) +
		geom_point() +
		geom_line() +
	  labs(x = "Iterations", y = "Mean Absolute Error") +
	  scale_x_continuous(breaks = pretty_breaks(n = length(Model$mae)))
	return(g)
}
