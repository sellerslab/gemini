 #' intialize_x
#' @description Initialize all x values to the value of \code{concordance}
#'
#' @param Model a Model object of class gemini.model
#' @param Input an Input object of class gemini
#' @param LFC an object within \code{Input} containing log-fold change values
#' @param guide.pair.annot an object within \code{Input} mapping guide pairs to individual genes.  See details for format.
#' @param concordance a numeric value to initialize x
#' @param pattern_split a character to strsplit the guide combinations found in \code{guide.pair.annot}
#' @param nc_gene a character naming the gene to use as a negative control
#' @param cores a numeric indicating the number of cores to use.  See \code{\link[gemini]{gemini_parallel}} default 1.
#' @param verbose default FALSE
#'
#' @details
#' guide.pair.annot must be (at least) of the following format:
#' \itemize{
#'  \item 1st column contains guide pairs, joined by \code{pattern_split}.
#'  \item 2nd column contains the name of gene mapping to the first guide in the 1st column.
#'  \item 3rd column contains the name of gene mapping to the second guide in the 1st column.
#' }
#'
#' @return a Model object of class gemini.model including new slots for x values
#'
#' @export
initialize_x <- function(Model,
						 Input, 
						 LFC = "LFC", 
						 guide.pair.annot = "guide.pair.annot",
						 concordance = 1, 
						 pattern_split = ";", 
						 nc_gene, 
						 cores = 1, 
						 verbose = F){
	# Check input
	stopifnot("gemini.model" %in% class(Model))
	
	# User message
	if(verbose) message("Initializing x")
	
	LFC <- Input[[LFC]]
	guide2gene <- Input[[guide.pair.annot]]

	# single guide to paired genes hashmap
	hashes = Sguide2Pguides_hash(guide2gene, pattern_split, cores = cores)

	# paired guides corresponding to nc_gene;nc_gene
	remove_seq1 = guide2gene[(guide2gene[,2] %in% nc_gene), 1]
	remove_seq2 = guide2gene[(guide2gene[,3] %in% nc_gene), 1]

	# guides designed for nc_gene
	paired_guide_remove_seq1 = as.character(hashes$paired_guide[remove_seq1,1]) %>%
		unique()
	paired_guide_remove_seq2 = as.character(hashes$paired_guide[remove_seq2,2]) %>%
		unique()
	paired_guide_remove = c(paired_guide_remove_seq1,paired_guide_remove_seq2)

	# remove nc_gene guides from hashes
	hashes$hash1 = hashes$hash1[!is.element(names(hashes$hash1),paired_guide_remove)]
	hashes$hash2 = hashes$hash2[!is.element(names(hashes$hash2),paired_guide_remove)]
	Model$hashes_x <- hashes

	# single guide concordance
	if(length(concordance)==1){
	  Model$x_seq1 = x_seq1 = rep(concordance,length(names(hashes$hash1))) %>% set_names(names(hashes$hash1))
	  Model$x_seq2 = x_seq2 = rep(concordance,length(names(hashes$hash2))) %>% set_names(names(hashes$hash2))
	  Model$x2_seq1 = x_seq1^2
	  Model$x2_seq2 = x_seq2^2
	}else if(length(concordance)==2 & is.list(concordance)){
	  Model$x_seq1 = x_seq1 = concordance[[1]][names(hashes$hash1)] %>% as.numeric() %>% set_names(names(hashes$hash1))
	  Model$x_seq2 = x_seq2 = concordance[[2]][names(hashes$hash2)] %>% as.numeric() %>% set_names(names(hashes$hash2))
	  Model$x2_seq1 = x_seq1^2
	  Model$x2_seq2 = x_seq2^2
	}


	# hash
	# genes corresponding to x_seq1
	x_seq1_genes = lapply(X = names(x_seq1), FUN = function(x) {
		gihj = hashes$hash1[[x]][1]
		g = guide2gene[match(gihj, guide2gene[,1]),2]
		return(g)
	}) %>% unlist() %>% set_names(names(x_seq1))

	# genes corresponding to x_seq2
	x_seq2_genes = lapply(names(x_seq2), function(x) {
		gihj = hashes$hash2[[x]][1]
		g = guide2gene[match(gihj, guide2gene[,1]),3]
		return(g)
	}) %>% unlist() %>% set_names(names(x_seq2))

	# paied guides concordance
	xx_names = guide2gene[!(guide2gene[,2] %in% nc_gene | guide2gene[,3] %in% nc_gene),1]
	Model$xx = xx = rep(concordance,length(xx_names)) %>% set_names(xx_names)
	Model$xx2 = xx^2
	Model$hashes_x$gene_hash1 <- x_seq1_genes
	Model$hashes_x$gene_hash2 <- x_seq2_genes


	# output
	return(Model)
}
