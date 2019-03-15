#' update_mae
#'
#' @description Calculate mean absolute error across all guide combinations.
#'
#' @param Model a Model object of class gemini.model
#' @param Input an Input object of class gemini.input
#' @param LFC an object within \code{Input} containing log-fold change values
#' @param guide.pair.annot an object within \code{Input} mapping guide pairs to individual genes.
#' @param nc_gene a character naming the gene to use as a negative control
#' @param pattern_join a character to join the gene combinations found in \code{guide.pair.annot}. default ';'
#' @param verbose default FALSE
#'
#' @return An object of class gemini.model
#'
#' @export
update_mae  <- function(Model, 
                        Input, 
                        LFC = "LFC", 
                        guide.pair.annot = "guide.pair.annot", 
                        nc_gene, 
                        pattern_join = ";",
                        verbose = F){
	stopifnot("gemini.model" %in% class(Model))
	
	# Calculate MAE
  if(verbose) message("Updating mae")
	LFC = Input[[LFC]]
	guide2gene = Input[[guide.pair.annot]]

	# Find each level of information - guide pair/guide, gene pair/gene
	gih = rownames(LFC)
	gi = Model$hashes_x$paired_guide[gih,1]
	hj = Model$hashes_x$paired_guide[gih,2]
	g = guide2gene[match(gih,guide2gene[,1]),2]
	h = guide2gene[match(gih,guide2gene[,1]),3]
	gh = apply(cbind(g,h), 1, function(x) paste(sort(x),collapse = pattern_join)) 
	
	# Initialize empty error matrix
	error = matrix(0, nrow = nrow(LFC), ncol = ncol(LFC), dimnames = list(rownames(LFC), colnames(LFC)))

	# Calculate error for each non-negative control pair
	ind = (!g %in% nc_gene) & (!h %in% nc_gene)
	error[ind,] = LFC[ind,] - Model$x_seq1[gi[ind]]*Model$y[g[ind],] -
		Model$x_seq2[hj[ind]]*Model$y[h[ind],] - Model$xx[gih[ind]]*Model$s[gh[ind],]
	
	# Calculate error for each pairs with 1 NC
	ind = (!g %in% nc_gene) & (h %in% nc_gene)
	error[ind,] = LFC[ind,] - Model$x_seq1[gi[ind]]*Model$y[g[ind],]

	ind = (g %in% nc_gene) & (!h %in% nc_gene)
	error[ind,] = LFC[ind,] - Model$x_seq2[hj[ind]]*Model$y[h[ind],]

	# Do not calculate error for pairs with 2 NCs
	ind = (g %in% nc_gene & h %in% nc_gene)
	error[ind,] = NA

	# output
	Model$mae = c(Model$mae, mean(abs(error), na.rm = T))
	Model$residuals = error
	return(Model)
}
