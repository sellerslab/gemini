#' update_y_pb
#' @description Update values of y using data from \code{Input} and current values of other parameters.
#'
#' @param Model a Model object of class gemini.model
#' @param Input an Input object of class gemini.input
#' @param LFC an object within \code{Input} containing log-fold change values
#' @param guide.pair.annot an object within \code{Input} mapping guide pairs to individual genes.
#' @param mean_y numeric indicating prior mean of y
#' @param sd_y numeric indicating prior sd of y
#' @param nc_gene a character naming the gene to use as a negative control
#' @param pattern_join a character to join the gene combinations found in \code{guide.pair.annot}
#' @param verbose default FALSE
#'
#' @return An object of class gemini.model
#'
#' @import pbmcapply
#' @export
update_y_pb <- function(Model,
						Input,
						LFC = "LFC",
						guide.pair.annot = "guide.pair.annot",
						mean_y = 0,
						sd_y = 10,
						nc_gene,
						pattern_join = ';',
						verbose = F) {
	if (verbose) {
		message("Updating y...")
		message("\tUsing 1 core (serial).")
		tstart = Sys.time()
		pb <- progressBar(min = 1, max = length(rownames(Model$y)))
	}
	
	LFC <- Input[[LFC]]
	guide2gene <- Input[[guide.pair.annot]]
	
	# genes corresponding to x_seq1
	x_seq1_genes = Model$hashes_x$gene_hash1
	
	# genes corresponding to x_seq2
	x_seq2_genes = Model$hashes_x$gene_hash2
	
	# mean of gamma distribution
	tau = Model$alpha / Model$beta
	
	# update y
	for (i in rownames(Model$y)) {
		gi_seq1 = names(x_seq1_genes)[x_seq1_genes == i]
		gi_seq2 = names(x_seq2_genes)[x_seq2_genes == i]
		gi_seq1_hj = lapply(gi_seq1, function(x)
			unlist(Model$hashes_x$hash1[x])) %>% set_names(gi_seq1)
		gi_seq2_hj = lapply(gi_seq2, function(x)
			unlist(Model$hashes_x$hash2[x])) %>% set_names(gi_seq2)
		
		# calculating updates for each x_seq1
		numerator = numeric(ncol(LFC)) %>% setNames(colnames(LFC))
		denominator = numeric(ncol(LFC)) %>% setNames(colnames(LFC))
		for (j in gi_seq1) {
			ind_nc = guide2gene[match(gi_seq1_hj[[j]], guide2gene[, 1]), 3] %in% nc_gene
			hj = Model$hashes_x$paired_guide[gi_seq1_hj[[j]], 2][!ind_nc]
			h = x_seq2_genes[hj]
			gh = apply(cbind(rep(i, length(h)), h), 1, function(x)
				paste(sort(x), collapse = pattern_join)) # gene g paired with genes h
			# numerator for each x_seq1
			numerator_seq1 = LFC[gi_seq1_hj[[j]][!ind_nc], ] - Model$x_seq2[hj] *
				Model$y[h, ] - Model$xx[gi_seq1_hj[[j]][!ind_nc]] * Model$s[gh, ]
			numerator_seq1 = (Model$x_seq1[j] * tau[gi_seq1_hj[[j]][!ind_nc], ]) *
				numerator_seq1
			numerator_seq1 = colSums(as.matrix(numerator_seq1), na.rm = T)
			numerator_seq1_nc = (Model$x_seq1[j] * tau[gi_seq1_hj[[j]][ind_nc], ]) *
				LFC[gi_seq1_hj[[j]][ind_nc], ]
			numerator_seq1_nc = colSums(as.matrix(numerator_seq1_nc), na.rm = T)
			numerator_seq1 = numerator_seq1 + numerator_seq1_nc
			numerator = numerator_seq1 + numerator
			# denominator for each x_seq1
			denominator_seq1 = colSums(as.matrix(Model$x2_seq1[j] * tau[gi_seq1_hj[[j]], ]), na.rm = T)
			denominator = denominator_seq1 + denominator
		}
		
		# calculating updates for each x_seq2
		for (j in gi_seq2) {
			ind_nc = guide2gene[match(gi_seq2_hj[[j]], guide2gene[, 1]), 2] %in% nc_gene
			hj = Model$hashes_x$paired_guide[gi_seq2_hj[[j]], 1][!ind_nc]
			h = x_seq1_genes[hj]
			gh = apply(cbind(rep(i, length(h)), h), 1, function(x)
				paste(sort(x), collapse = pattern_join)) # gene g paired with genes h
			# numerator for each x_seq2
			numerator_seq2 = LFC[gi_seq2_hj[[j]][!ind_nc], ] - Model$x_seq1[hj] *
				Model$y[h, ] - Model$xx[gi_seq2_hj[[j]][!ind_nc]] * Model$s[gh, ]
			numerator_seq2 = (Model$x_seq2[j] * tau[gi_seq2_hj[[j]][!ind_nc], ]) *
				numerator_seq2
			numerator_seq2 = colSums(as.matrix(numerator_seq2), na.rm = T)
			numerator_seq2_nc = (Model$x_seq2[j] * tau[gi_seq2_hj[[j]][ind_nc], ]) *
				LFC[gi_seq2_hj[[j]][ind_nc], ]
			numerator_seq2_nc = colSums(as.matrix(numerator_seq2_nc), na.rm = T)
			numerator_seq2 = numerator_seq2 + numerator_seq2_nc
			numerator = numerator_seq2 + numerator
			# denominator for each x_seq2
			denominator_seq2 = colSums(as.matrix(Model$x2_seq2[j] * tau[gi_seq2_hj[[j]], ]), na.rm = T)
			denominator = denominator_seq2 + denominator
		}
		
		# update progress bar
		if (verbose)
			setTxtProgressBar(pb = pb, value = getTxtProgressBar(pb) + 1)
		
		# adding prior values to the numerator and denominator
		numerator = mean_y / (sd_y ^ 2) + numerator
		denominator = 1 / (sd_y ^ 2) + denominator
		
		# updating y and y2
		Model$y[i, ] = numerator / denominator
		Model$y2[i, ] = Model$y[i, ] ^ 2 + 1 / denominator
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
