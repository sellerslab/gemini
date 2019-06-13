#' update_mae
#'
#' @description Calculate mean absolute error across all guide combinations.
#'
#' @param Model a Model object of class gemini.model
#' @param verbose default FALSE
#'
#' @return An object of class gemini.model
#' 
#' @examples 
#' Model <- gemini::Model
#' Model %<>% update_mae()
#'
#' @export
update_mae  <- function(Model, 
                        verbose = FALSE){
	stopifnot("gemini.model" %in% class(Model))
	
	# Calculate MAE
	if(verbose) message("Updating mae")
	
	Input <- Model$Input
	LFC = Input[[Model$LFC.name]]
	guide2gene = Input[[Model$guide.pair.annot]]

	# Find each level of information - guide pair/guide, gene pair/gene
	gihj = rownames(LFC)
	gi = Model$hashes_x$paired_guide[gihj,1]
	hj = Model$hashes_x$paired_guide[gihj,2]
	g = guide2gene[match(gihj,guide2gene[,1]),2]
	h = guide2gene[match(gihj,guide2gene[,1]),3]
	gh = apply(cbind(g,h), 1, function(x) paste(sort(x),collapse = Model$pattern_join)) 
	
	# Initialize empty error matrix
	error = matrix(0, nrow = nrow(LFC), ncol = ncol(LFC), dimnames = list(rownames(LFC), colnames(LFC)))

	# Calculate error for each non-negative control pair
	ind = (!g %in% Model$nc_gene) & (!h %in% Model$nc_gene)
	error[ind,] = LFC[ind,] - Model$x[gi[ind]]*Model$y[g[ind],] -
		Model$x[hj[ind]]*Model$y[h[ind],] - Model$xx[gihj[ind]]*Model$s[gh[ind],]
	
	# Calculate error for each pairs with 1 NC
	ind = (!g %in% Model$nc_gene) & (h %in% Model$nc_gene)
	error[ind,] = LFC[ind,] - Model$x[gi[ind]]*Model$y[g[ind],]

	ind = (g %in% Model$nc_gene) & (!h %in% Model$nc_gene)
	error[ind,] = LFC[ind,] - Model$x[hj[ind]]*Model$y[h[ind],]

	# Do not calculate error for pairs with 2 NCs
	ind = (g %in% Model$nc_gene & h %in% Model$nc_gene)
	error[ind,] = NA

	# output
	Model$mae = c(Model$mae, mean(abs(error), na.rm = TRUE))
	Model$residuals = error
	return(Model)
}
