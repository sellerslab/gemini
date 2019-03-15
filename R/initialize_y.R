#' initialize_y
#' @description Initialize all y values from guide pairs including a negative control.
#'
#' @param Model an object of class gemini.model
#' @param Input an object of class gemini.input
#' @param LFC an object within \code{Input} containing log-fold change values
#' @param guide.pair.annot an object within \code{Input} containing guide to gene annotations
#' @param nc_gene a character naming the gene to use as a negative control
#' @param verbose default FALSE
#' @param cores a numeric indicating the number of cores to use.  See details in gemini_parallel. (default=1)
#'
#' @return a Model object of class gemini.model including new slots for y values
#'
#' @export
initialize_y <- function(Model,
						 Input, 
						 LFC = "LFC", 
						 guide.pair.annot = "guide.pair.annot",
						 nc_gene, 
						 verbose = F,
						 cores = 1){
	# Check input
	stopifnot("gemini.model" %in% class(Model))
	if(verbose) message("Initializing y")
	
	# Create inputs
	guide2gene <- Input[[guide.pair.annot]]
	LFC <- Input[[LFC]]

	# gene pairs containing nc_gene
	guide2nc_gene = guide2gene[guide2gene[,2] %in% nc_gene | guide2gene[,3] %in% nc_gene,]

	# genes mapped to guides
	hash_nc = Sgene2Pguides_hash(guide2nc_gene, cores)
	hash_others = Sgene2Pguides_hash(guide2gene, cores)

	# gene list
	genes = unique(c(guide2gene[,2], guide2gene[,3]))

	# removing nc_gene from gene list
	genes = genes[!is.element(genes,nc_gene)]

	# define matrix y
	y = matrix(nrow = length(genes), ncol = ncol(LFC), dimnames = list(genes, colnames(LFC)))

	# estimate y
	for (i in genes){
		ind = hash_nc[[i]]
		if(length(ind)>0){
			y[i,] = apply(as.matrix(LFC[ind,]), 2, function(x) median(x, na.rm = T))
		}else{
			ind = hash_others
			y[,i] = apply(as.matrix(LFC[ind,]), 2, function(x) median(x, na.rm = T))
		}
	}

	# output
	Model$y <- y
	Model$y2 <- y^2

	return(Model)
	#return(list(y=y, y2=y^2, hash_y = Sgene2Pguides_hash(guide2gene)[genes])) # no need for hash_y
}
