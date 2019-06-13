#' update_y_pb
#' @description Update values of y using data from \code{Input} and current values of other parameters.
#'
#' @param Model a Model object of class gemini.model
#' @param mean_y numeric indicating prior mean of y
#' @param sd_y numeric indicating prior sd of y
#' @param verbose default FALSE
#'
#' @return An object of class gemini.model
#'
#' @import pbmcapply
#' @importFrom stats setNames
#' @importFrom utils getTxtProgressBar
#' @importFrom utils setTxtProgressBar
#' 
#' @examples
#' Model <- gemini::Model
#' Model %<>% update_y_pb()
#' 
#' @export
update_y_pb <- function(Model,
						mean_y = 0,
						sd_y = 10,
						verbose = FALSE) {
	if (verbose) {
		message("Updating y...")
		message("\tUsing 1 core (serial).")
		tstart = Sys.time()
		pb <- progressBar(min = 1, max = length(rownames(Model$y)))
	}
	
	Input <- Model$Input
	LFC <- Input[[Model$LFC.name]]
	guide2gene <- Input[[Model$guide.pair.annot]]
	
	# genes corresponding to x_seq1
	x_seq_genes = Model$hashes_x$gene_hash

	# mean of gamma distribution
	tau = Model$alpha / Model$beta
	
	# update y
	for (g in rownames(Model$y)) {
		gi_seq = names(x_seq_genes)[x_seq_genes == g]
		gi_seq_hj = lapply(gi_seq, function(x)
			unlist(Model$hashes_x$hash[x])) %>% set_names(gi_seq)
		
		# calculating updates for each x_seq
		numerator = denominator = numeric(ncol(LFC)) %>% setNames(colnames(LFC))
		
		for (i in gi_seq) {
			# identify which position the partner gene (h) is in
			h.col <- vapply(gi_seq_hj[[i]], function(s){
				ss = strsplit(s, split = Model$pattern_split, fixed = TRUE)[[1]]
				return(which(ss != i))
			}, numeric(1))
			
			# identify cases in which the partner gene is nc_gene
			ind_nc = guide2gene[match(gi_seq_hj[[i]], guide2gene[, 1]), 3]
			ind_nc[h.col==1] = guide2gene[match(gi_seq_hj[[i]], guide2gene[, 1]), 2][h.col==1]
			ind_nc = ind_nc %in% Model$nc_gene
			
			# identify partner guides, filtering out nc_gene guides
			hj = Model$hashes_x$paired_guide[gi_seq_hj[[i]], 2]
			hj[h.col==1] = Model$hashes_x$paired_guide[gi_seq_hj[[i]], 1][h.col==1]
			hj = hj[!ind_nc]
			
			# gene g paired with genes h
			h = x_seq_genes[hj]
			gh = apply(cbind(rep(g, length(h)), h), 1, function(x)
				paste(sort(x), collapse = Model$pattern_join)) 
			
			# numerator for each guide i
			numerator_i = LFC[gi_seq_hj[[i]][!ind_nc], ] - Model$x[hj] *
				Model$y[h, ] - Model$xx[gi_seq_hj[[i]][!ind_nc]] * Model$s[gh, ]
			numerator_i = (Model$x[i] * tau[gi_seq_hj[[i]][!ind_nc], ]) *
				numerator_i
			if (sum(!ind_nc)>1 | sum(!ind_nc)==0){
				numerator_i = colSums(as.matrix(numerator_i), na.rm = TRUE)
			} else{
				numerator_i = colSums(as.matrix(t(numerator_i)), na.rm = TRUE)
			}
			numerator_nc_i = (Model$x[i] * tau[gi_seq_hj[[i]][ind_nc], ]) *
				LFC[gi_seq_hj[[i]][ind_nc], ]
			if (sum(ind_nc)>1 | sum(ind_nc)==0){
				numerator_nc_i = colSums(as.matrix(numerator_nc_i), na.rm = TRUE)	
			} else{
				numerator_nc_i = colSums(as.matrix(t(numerator_nc_i)), na.rm = TRUE)	
			}
			numerator_i = numerator_i + numerator_nc_i
			numerator = numerator_i + numerator
			
			# denominator for each guide i
			if (length(gi_seq_hj[[i]])>1){
				denominator_i = colSums(as.matrix(Model$x2[i] * tau[gi_seq_hj[[i]], ]), na.rm = TRUE)	
			} else{
				denominator_i = colSums(as.matrix(t(Model$x2[i] * tau[gi_seq_hj[[i]], ])), na.rm = TRUE)
			}
			#denominator_i = colSums(as.matrix(Model$x2[i] * tau[gi_seq_hj[[i]], ]), na.rm = TRUE)
			denominator = denominator_i + denominator
		}
		
		# update progress bar
		if (verbose)
			setTxtProgressBar(pb = pb, value = getTxtProgressBar(pb) + 1)
		
		# adding prior values to the numerator and denominator
		numerator = mean_y / (sd_y ^ 2) + numerator
		denominator = 1 / (sd_y ^ 2) + denominator
		
		# updating y and y2
		Model$y[g, ] = numerator / denominator
		Model$y2[g, ] = Model$y[g, ] ^ 2 + 1 / denominator
	}
	
	
	# output
	if (verbose) {
		close(pb)
		tend = Sys.time()
		tdiff = difftime(tend, tstart)
		message("\tCompleted update of y.")
		message("\tTime to completion: ",
				round(tdiff, digits = 3),
				' ',
				units(tdiff))
	}
	return(Model)
}
