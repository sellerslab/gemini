#' update_tau_pb
#' @description Update parameters of tau using data from \code{Input} and current values of other parameters.
#'
#' @param Model a Model object of class gemini.model
#' @param cores a numeric indicating the number of cores to use.  See \code{\link[gemini]{gemini_parallelization}} for details.  (default=1).
#' @param verbose default FALSE
#'
#' @return An object of class gemini.model
#'
#' @import pbmcapply
#' @import magrittr
#'
#' @examples 
#' Model <- gemini::Model
#' Model %<>% update_tau_pb()
#'
#' @export
update_tau_pb <- function(Model,
						  cores = 1,
						  verbose = FALSE) {
	if (verbose) {
		message("Updating tau...")
		message("\tUsing ", cores, ' core(s).')
		tstart = Sys.time()
	}
	
	Input <- Model$Input
	LFC <- Input[[Model$LFC.name]]
	guide2gene <- Input[[Model$guide.pair.annot]]
	
	# mean of gamma distribution
	tau = Model$alpha / Model$beta
	
	# define loop function
	tau_loop <-
		function(gihj,
				 Model,
				 guide2gene,
				 LFC) {
			gi = Model$hashes_x$paired_guide[gihj, 1]
			hj = Model$hashes_x$paired_guide[gihj, 2]
			g = guide2gene[match(gihj, guide2gene[, 1]), 2]
			h = guide2gene[match(gihj, guide2gene[, 1]), 3]
			gh = paste(sort(c(g, h)), collapse = Model$pattern_join)
			
			# calculating beta*
			gnc = g %in% Model$nc_gene
			hnc = h %in% Model$nc_gene
			if (!(gnc) & !(hnc)) {
				beta_star = LFC[gihj,] ^ 2 - 
					2 * LFC[gihj,] * (Model$x[gi] * Model$y[g,] + Model$x[hj] * Model$y[h,] + Model$xx[gihj] * Model$s[gh,]) + 
					Model$x2[gi] * Model$y2[g,] + 
					2 * Model$x[gi] * Model$y[g,] * (Model$x[hj] * Model$y[h,] + Model$xx[gihj] * Model$s[gh,]) + 
					Model$x2[hj] * Model$y2[h,] + Model$xx2[gihj] * Model$s2[gh,] + 
					2 * Model$x[hj] * Model$y[h,] * Model$xx[gihj] * Model$s[gh,]
			} else if (gnc & !(hnc)) {
				beta_star = LFC[gihj,] ^ 2 + 
					Model$x2[hj] * Model$y2[h,] - 
					2 * LFC[gihj,] * Model$x[hj] * Model$y[h,]
			} else if (!(gnc) & hnc) {
				beta_star = LFC[gihj,] ^ 2 + 
					Model$x2[gi] * Model$y2[g,] - 
					2 * LFC[gihj,] * Model$x[gi] * Model$y[g,]
			} else {
				beta_star = LFC[gihj,] ^ 2
			}
			
			# updating alpha and beta
			alpha_gihj = Model$prior_shape + 0.5
			beta_gihj = Model$beta_prior[gihj,] + 0.5 * beta_star
			return(list(alpha = alpha_gihj, beta = beta_gihj))
		}
	
	if(verbose){
	    res <-
	        pbmcapply::pbmclapply(
	            X = rownames(Model$alpha),
	            FUN = tau_loop,
	            Model = Model,
	            guide2gene = guide2gene,
	            LFC = LFC,
	            mc.cores = cores
	        )
	}else{
	    res <-
	        parallel::mclapply(
	            X = rownames(Model$alpha),
	            FUN = tau_loop,
	            Model = Model,
	            guide2gene = guide2gene,
	            LFC = LFC,
	            mc.cores = cores
	        )
	}
	
	
	Model$beta[,] <- lapply(res, magrittr::extract, "beta") %>%
		unlist(recursive = FALSE, use.names = FALSE) %>%
		do.call(rbind, .)
	
	Model$alpha[,] <- lapply(res, magrittr::extract, "alpha") %>%
		unlist(recursive = FALSE, use.names = FALSE) %>%
		do.call(rbind, .)
	
	# output
	if (verbose) {
		tend = Sys.time()
		tdiff = difftime(tend, tstart)
		message("\tCompleted update of tau.")
		message("\tTime to completion: ",
				round(tdiff, digits = 3),
				' ',
				units(tdiff))
	}
	return(Model)
}
