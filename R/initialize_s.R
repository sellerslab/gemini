#' initialize_s
#' @description Initialize s values using initialized y values and data from \code{Input}
#'
#' @param Model an object of class gemini.model
#' @param Input an object of class gemini.input
#' @param LFC an object within \code{Input} containing log-fold change values
#' @param guide.pair.annot an object within \code{Input} containing guide to gene annotations
#' @param pattern_join a character to join the gene combinations found in \code{guide.pair.annot}
#' @param nc_gene a character naming the gene to use as a negative control
#' @param cores a numeric indicating the number of cores to use.  See \code{\link[gemini]{gemini_parallel}} default 1.
#' @param verbose default FALSE
#'
#' @return a Model object of class gemini.model including new slots for s values
#'
#' @export
initialize_s <- function(Model, 
						 Input, 
						 LFC, 
						 guide.pair.annot = "guide.pair.annot",
						 nc_gene,
						 pattern_join = ";",
						 cores = 1,
						 verbose = F){
	# Check input
	stopifnot("gemini.model" %in% class(Model))
	
	# User message
	if(verbose) message("Initializing s")
	
	guide2gene <- Input[[guide.pair.annot]]
	LFC <- Input[[LFC]]

	# calculate initial y
	y = Model$y

	# remove control gene
	guide2gene = guide2gene[!(guide2gene[,2] %in% nc_gene | guide2gene[,3] %in% nc_gene),]

	# find unique gene pairs
	guide2gene$Key <- apply(guide2gene[,c(2,3)], 1, function(x) paste(sort(x),collapse = pattern_join))
	Pairs = unique(guide2gene$Key)

	# paired genes mapped to paired guides
	hash <- parallel::mclapply(Pairs, function(x){
		as.character(guide2gene[guide2gene$Key == x,1])
	}, mc.cores = cores) %>% set_names(Pairs)
	
	# Initialize s using initialized y values
	res = parallel::mclapply(Pairs, function(i){
		ind = unlist(hash[i])
		vec = apply(as.matrix(LFC[ind,]), 2, function(x) median(x, na.rm = T))
		g = strsplit(i,pattern_join)[[1]][1]
		h = strsplit(i,pattern_join)[[1]][2]
		vec = vec - y[g,] - y[h,]
		return(vec)
	}, mc.cores = cores)

	s = matrix(unlist(res), byrow = T, nrow = length(Pairs), ncol = ncol(LFC), dimnames = list(Pairs,colnames(LFC)))
	
	# save output
	Model$s <- s
	Model$s2 <- s^2
	Model$hash_s <- hash

	return(Model)
}
