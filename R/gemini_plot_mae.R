#' MAE plot
#' @name gemini_plot_mae
#'
#' @description Plots the average MAE of gemini model at each iteration step.
#'
#' @param Model a gemini.model object
#'
#' @import ggplot2
#' @importFrom scales pretty_breaks
#' @export
gemini_plot_mae <- function(Model){
	stopifnot("gemini.model" %in% class(Model))
	stopifnot("mae" %in% names(Model))

	g = ggplot(data = data.frame(iter = (1:length(Model$mae))-1, mae = Model$mae),
					aes(x = iter, y = mae)) +
		geom_point() +
		geom_line() +
	  labs(x = "Iterations", y = "Mean Absolute Error") +
	  scale_x_continuous(breaks = pretty_breaks(n = length(Model$mae)))
	return(g)
}
