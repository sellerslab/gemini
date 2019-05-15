#' intialize_x
#' @description Initialize all x values to the value of \code{concordance}
#'
#' @param Model a Model object of class gemini.model
#' @param concordance a numeric value to initialize x
#' @param cores a numeric indicating the number of cores to use.  See \code{\link[gemini]{gemini_parallel}} default 1.
#' @param verbose default FALSE
#'
#' @note As there is much hashing involved in this function, this tends to be computationally intensive.
#' As such, we have enabled parallelization of most hash steps, but this may still be rate-limited by the 
#' amount of memory consumed. 
#' 
#' @return a Model object of class gemini.model including new slots for x values and internal-use hashes
#'
#' @export
initialize_x <- function(Model,
						 concordance = 1,
						 cores = 1,
						 verbose = FALSE) {
	# Check input
	stopifnot("gemini.model" %in% class(Model))
	
	# User message
	if (verbose)
		message("Initializing x")
	
	Input <- Model$Input
	LFC <- Input[[Model$LFC.name]]
	guide2gene <- Input[[Model$guide.pair.annot]]
	
	# single guide to paired genes hashmap
	hash = Sguide2Pguides_hash(guide2gene, Model$pattern_split, cores = cores)
	
	# paired guides corresponding to nc_gene;nc_gene
	remove_seq1 = guide2gene[(guide2gene[, 2] %in% Model$nc_gene), 1]
	remove_seq2 = guide2gene[(guide2gene[, 3] %in% Model$nc_gene), 1]
	
	# guides designed for nc_gene
	paired_guide_remove_seq1 = as.character(hash$paired_guide[remove_seq1, 1]) %>%
		unique()
	paired_guide_remove_seq2 = as.character(hash$paired_guide[remove_seq2, 2]) %>%
		unique()
	paired_guide_remove = c(paired_guide_remove_seq1, paired_guide_remove_seq2)
	
	# remove nc_gene guides from hashes
	hash$hash = hash$hash[!is.element(names(hash$hash), paired_guide_remove)]
	Model$hashes_x <- hash
	
	# single guide concordance
	if (length(concordance) == 1 & is.numeric(concordance)) {
		Model$x = x_seq = rep(concordance, length(names(hash$hash))) %>% set_names(names(hash$hash))
		Model$x2 = x_seq ^ 2
	} else{
		stop(
			"No single concordance value specified. Please specify a single concordance value for x."
		)
	}
	
	# hash
	# genes corresponding to x_seq
	x_seq_genes = parallel::mclapply(
		X = names(x_seq),
		FUN = function(x) {
			gihj = hash$hash[[x]][1]
			ij = strsplit(gihj, split = Model$pattern_split, fixed = TRUE)[[1]]
			gene.col = which(ij == x) + 1
			g = guide2gene[match(gihj, guide2gene[, 1]), gene.col]
			return(g)
		},
		mc.cores = cores
	) %>% unlist() %>% set_names(names(x_seq))
	
	Model$hashes_x$gene_hash <- x_seq_genes

	# paired guides concordance
	xx_names = guide2gene[!(guide2gene[, 2] %in% Model$nc_gene |
								guide2gene[, 3] %in% Model$nc_gene), 1]
	Model$xx = xx = rep(concordance, length(xx_names)) %>% set_names(xx_names)
	Model$xx2 = xx ^ 2
	
	# output
	return(Model)
}
