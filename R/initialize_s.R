#' initialize_s
#' @description Initialize s values using initialized y values and data from \code{Input}
#'
#' @param Model an object of class gemini.model
#' @param cores a numeric indicating the number of cores to use.  See \code{\link[gemini]{gemini_parallelization}} default 1.
#' @param verbose default FALSE
#'
#' @return a Model object of class gemini.model including new slots for s values
#'
#' @export
#' 
#' @examples
#' data("Model", package = "gemini")
#' Model <- initialize_s(Model)
initialize_s <- function(Model,
						 cores = 1,
						 verbose = FALSE) {
	# Check input
	stopifnot("gemini.model" %in% class(Model))
	
	# User message
	if (verbose)
		message("Initializing s")
	
	Input <- Model$Input
	guide2gene <- Input[[Model$guide.pair.annot]]
	LFC <- Input[[Model$LFC.name]]
	
	# calculate initial y
	y = Model$y
	
	# remove control gene
	guide2gene = guide2gene[!(guide2gene[, 2] %in% Model$nc_gene |
							  	guide2gene[, 3] %in% Model$nc_gene), ]
	
	# find unique gene pairs
	guide2gene$Key <-
		apply(guide2gene[, c(2, 3)], 1, function(x)
			paste(sort(x), collapse = Model$pattern_join))
	Pairs = unique(guide2gene$Key)
	
	# paired genes mapped to paired guides
	hash <- parallel::mclapply(Pairs, function(x) {
		as.character(guide2gene[guide2gene$Key == x, 1])
	}, mc.cores = cores) %>% set_names(Pairs)
	
	# Initialize s using initialized y values
	res = parallel::mclapply(Pairs, function(i) {
		ind = unlist(hash[i])
		if (length(ind)>1){
			vec = apply(as.matrix(LFC[ind, ]), 2, function(x)
				median(x, na.rm = TRUE))
		} else{
			vec = apply(as.matrix(t(LFC[ind, ])), 2, function(x)
				median(x, na.rm = TRUE))	
		}
		g = strsplit(i, Model$pattern_join)[[1]][1]
		h = strsplit(i, Model$pattern_join)[[1]][2]
		vec = vec - y[g, ] - y[h, ]
		return(vec)
	}, mc.cores = cores)
	
	s = matrix(
		unlist(res),
		byrow = TRUE,
		nrow = length(Pairs),
		ncol = ncol(LFC),
		dimnames = list(Pairs, colnames(LFC))
	)
	
	# save output
	Model$s <- s
	Model$s2 <- s ^ 2
	Model$hash_s <- hash
	
	return(Model)
}
