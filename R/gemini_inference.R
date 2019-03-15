#' gemini_inference
#'
#' @description
#' Estimate the posterior using a variational inference technique. Inference is performed
#' through an iterative process until convergence.
#'
#' @param Model an object of class gemini.model
#' @param mean_x a numeric indicating prior mean of x (default=1).
#' @param sd_x a numeric indicating prior sd of x (default=1).
#' @param mean_xx a numeric indicating prior mean of xx (default=1).
#' @param sd_xx a numeric indicating prior sd of xx (default=1)
#' @param mean_y a numeric indicating prior mean of y (default=0)
#' @param sd_y a numeric indicating prior sd of y (default=10).
#' @param mean_s a numeric indicating prior mean of s(default=0)
#' @param sd_s a numeric indicating prior sd of s (default=10)
#' @param n_iterations a numeric indicating the maximum number of iterations (default=20).
#' @param threshold a numeric indicating the threshold of change in MAE at which to stop the iterative process (default=0.001).
#' @param cores a numeric indicating the number of cores to use.  See details in gemini_parallel. (default=1)
#' @param force_results a logical indicating if the CAVI algorithm should be halted if non-convergence is detected. (default=FALSE)
#' @param verbose default FALSE
#' @param save_iterations for especially large libraries that require long computations,
#' saves the latest iteration of each update.  default FALSE
#'
#' @details
#' GEMINI uses the following parameters, which are described in Zamanighomi et al. and translated here for clarity:
#' \itemize{
#' \item \strong{y:} individual gene effect
#' \item \strong{s:} combination effect
#' \item \strong{x:} screen variation corresponding to individual guides
#' \item \strong{xx:} screen variation corresponding to paired guides
#' }
#'
#' Default parameters may need to be changed if convergence is not achieved. See README for
#' more details.
#'
#' @return a gemini.model object with estimated posteriors
#'
#' @export
gemini_inference <- function(Model,
							 n_iterations = 20,
							 mean_x = 1,
							 sd_x = 1,
							 mean_xx = 1,
							 sd_xx = 1,
							 mean_y = 0,
							 sd_y = 10,
							 mean_s = 0,
							 sd_s = 10,
							 threshold = 1e-3,
							 cores = 1,
							 force_results = F,
							 verbose = FALSE,
							 save_iterations = FALSE) {
	# perform n_iterations
	for (i in 1:n_iterations) {
		if (verbose)
			cat(paste("iteration:", i, "\n"))
		Model %<>%
			update_x_pb(
				Input = Model$Input,
				LFC = Model$LFC,
				guide.pair.annot = Model$guide.pair.annot,
				nc_gene = Model$nc_gene,
				mean_x = mean_x,
				sd_x = sd_x,
				mean_xx = mean_xx,
				sd_xx = sd_xx,
				pattern_join = Model$pattern_join,
				cores = cores,
				verbose = verbose
			) %>%
			update_y_pb(
				Input = Model$Input,
				LFC = Model$LFC,
				guide.pair.annot = Model$guide.pair.annot,
				nc_gene = Model$nc_gene,
				mean_y = mean_y,
				sd_y = sd_y,
				pattern_join = Model$pattern_join,
				verbose = verbose
			) %>%
			update_s_pb(
				Input = Model$Input,
				LFC = Model$LFC,
				guide.pair.annot = Model$guide.pair.annot,
				mean_s = mean_s,
				sd_s = sd_s,
				cores = cores,
				verbose = verbose
			) %>%
			update_tau_pb(
				Input = Model$Input,
				LFC = Model$LFC,
				guide.pair.annot = Model$guide.pair.annot,
				nc_gene = Model$nc_gene,
				pattern_join = Model$pattern_join,
				cores = cores,
				verbose = verbose
			) %>%
			update_mae(
				Input = Model$Input,
				LFC = Model$LFC,
				guide.pair.annot = Model$guide.pair.annot,
				nc_gene = Model$nc_gene,
				pattern_join = Model$pattern_join,
				verbose = verbose
			)
		
		# save the Model object after every iteration - useful for large screens with heavy computation
		if (!is.null(save_iterations)) {
			if (is.logical(save_iterations)) {
				save_flag = save_iterations
				save_file = "Model"
			} else if (is.character(save_iterations)) {
				save_flag = TRUE
				save_file = save_iterations
			}else{
				if (verbose)
					warning("Could not parse save_iterations argument.")
				save_flag = FALSE
			}
			if (isTRUE(save_flag)) {
				if (file.exists(paste0(save_file, "_", i - 1, ".rds"))) {
					file.remove(paste0(save_file, "_", i - 1, ".rds"))
				}
				saveRDS(Model, file = paste0(save_file, "_", i, ".rds"))
			}
		}
		
		# break loop if mae does not change more than the threshold
		if (verbose) {
			message("MAE: ", Model$mae[i + 1])
		}
		if (abs(Model$mae[i + 1] - Model$mae[i]) < threshold) {
			break
		}
		if (!force_results) {
			if (Model$mae[i + 1] > Model$mae[i]) {
				stop("Model not converging! Please adjust parameters and try again.")
			}
		}
	}
	# output
	return(Model)
}