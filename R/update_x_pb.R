#' update_x_p
#' @description Update values of x using data from \code{Input} and current values of other parameters.
#'
#' @param Model a Model object of class gemini.model
#' @param Input an Input object of class gemini.input
#' @param LFC an object within \code{Input} containing log-fold change values
#' @param guide.pair.annot an object within \code{Input} mapping guide pairs to individual genes.
#' @param mean_x a numeric indicating prior mean of x
#' @param sd_x a numeric indicating prior sd of x
#' @param mean_xx a numeric indicating prior mean of xx
#' @param sd_xx a numeric indicating prior sd of xx
#' @param nc_gene a character naming the gene to use as a negative control
#' @param pattern_join a character to join the gene combinations found in \code{guide.pair.annot}
#' @param cores a numeric indicating the number of cores to use.  See \code{\link{gemini_parallel}} for details.  (default=1).
#' @param verbose default FALSE
#'
#' @return An object of class gemini.model
#'
#' @import pbmcapply
#' @import magrittr
#'
#' @export
update_x_pb <- function(Model,
						Input,
						LFC = "LFC",
						guide.pair.annot = "guide.pair.annot",
						mean_x = 1,
						sd_x = 1,
						mean_xx = 1,
						sd_xx = 1,
						nc_gene,
						pattern_join = ";",
						cores = 1,
						verbose = F) {
	if (verbose) {
		message("Updating x...")
		message("\tUsing ", cores, " core(s).")
		tstart = Sys.time()
	}
	
	guide2gene <- Input[[guide.pair.annot]]
	LFC <- Input[[LFC]]
	
	# genes corresponding to x_seq1
	x_seq1_genes = Model$hashes_x$gene_hash1
	
	# genes corresponding to x_seq2
	x_seq2_genes = Model$hashes_x$gene_hash2
	
	# mean of gamma distribution
	tau = Model$alpha / Model$beta
	
	# Allow looping through sd_x if lists are provided
	if (length(sd_x) == 1) {
		sd_x1 = rep(sd_x, length(Model$hashes_x$hash1))
		sd_x2 = rep(sd_x, length(Model$hashes_x$hash2))
		names(sd_x1) <- names(Model$hashes_x$hash1)
		names(sd_x2) <- names(Model$hashes_x$hash2)
		sd_x = list(sd_x1 = sd_x1, sd_x2 = sd_x2)
	}
	if (length(mean_x) == 1) {
		mean_x1 = rep(mean_x, length(Model$hashes_x$hash1))
		mean_x2 = rep(mean_x, length(Model$hashes_x$hash2))
		names(mean_x1) <- names(Model$hashes_x$hash1)
		names(mean_x2) <- names(Model$hashes_x$hash2)
		mean_x = list(mean_x1 = mean_x1, mean_x2 = mean_x2)
	}
	if (length(mean_xx) == 1) {
		mean_xx = rep(mean_xx, nrow(Model$hashes_x$paired_guide))
		names(mean_xx) <- rownames(Model$hashes_x$paired_guide)
	}
	if (length(sd_xx) == 1) {
		sd_xx = rep(sd_xx, nrow(Model$hashes_x$paired_guide))
		names(sd_xx) <- rownames(Model$hashes_x$paired_guide)
	}
	
	# update x_seq1 and x2_seq1
	x1_loop <-
		function(i,
				 Model,
				 guide2gene,
				 LFC,
				 tau,
				 x_seq1_genes,
				 pattern_join,
				 nc_gene,
				 mean_x,
				 sd_x) {
			# get individual mean_xi and sd_xi for guide i
			mean_xi = mean_x[[1]][i]
			sd_xi = sd_x[[1]][i]
			# x_seq1 to gene names
			g = as.character(x_seq1_genes[i]) # gene g corresponding to guide i
			gihj = Model$hashes_x$hash1[[i]] # all guide pairs containing i
			h = guide2gene[match(gihj, guide2gene[, 1]), 3] # genes h
			hj = Model$hashes_x$paired_guide[match(gihj, rownames(Model$hashes_x$paired_guide)), 2] # guides for genes h
			gh = apply(cbind(rep(g, length(h)), h), 1, function(x)
				paste(sort(x), collapse = pattern_join)) # gene g paired with genes h
			# calculating numerator
			numerator = LFC[gihj[!h %in% nc_gene], ] - Model$x_seq2[hj[!h %in% nc_gene]] *
				Model$y[h[!h %in% nc_gene], ] - Model$xx[gihj[!h %in% nc_gene]] * Model$s[gh[!h %in% nc_gene], ]
			numerator = tau[gihj[!h %in% nc_gene], ] * numerator
			numerator = colSums(as.matrix(numerator), na.rm = T) +
				colSums(as.matrix(tau[gihj[h %in% nc_gene], ] * LFC[gihj[h %in% nc_gene], ]), na.rm = T)
			numerator = sum(Model$y[g, ] * numerator, na.rm = T)
			numerator = mean_xi / (sd_xi ^ 2) + numerator
			# calculating denominator
			denominator = sum(colSums(as.matrix(tau[gihj, ]), na.rm = T) * Model$y2[g, ], na.rm = T)
			denominator = 1 / (sd_xi ^ 2) + denominator
			x_seq = numerator / denominator
			x2_seq = x_seq ^ 2 + 1 / denominator
			return(list(x_seq = x_seq, x2_seq = x2_seq))
		}
	
	# update x_seq1 and x2_seq1
	x2_loop <-
		function(i,
				 Model,
				 guide2gene,
				 LFC,
				 tau,
				 x_seq2_genes,
				 pattern_join,
				 nc_gene,
				 mean_x,
				 sd_x) {
			# get individual mean_xi and sd_xi for guide i
			mean_xi = mean_x[[2]][i]
			sd_xi = sd_x[[2]][i]
			# x_seq2 to gene names
			g = as.character(x_seq2_genes[i]) # gene g
			gihj = Model$hashes_x$hash2[[i]]
			h = guide2gene[match(gihj, guide2gene[, 1]), 2] # genes h
			hj = Model$hashes_x$paired_guide[match(gihj, rownames(Model$hashes_x$paired_guide)), 1] # guides for genes h
			gh = apply(cbind(rep(g, length(h)), h), 1, function(x)
				paste(sort(x), collapse = pattern_join)) # gene g paired with genes h
			# calculating numerator
			numerator = LFC[gihj[!h %in% nc_gene], ] - Model$x_seq1[hj[!h %in% nc_gene]] *
				Model$y[h[!h %in% nc_gene], ] - Model$xx[gihj[!h %in% nc_gene]] * Model$s[gh[!h %in% nc_gene], ]
			numerator = tau[gihj[!h %in% nc_gene], ] * numerator
			numerator = colSums(as.matrix(numerator), na.rm = T) +
				colSums(as.matrix(tau[gihj[h %in% nc_gene], ] * LFC[gihj[h %in% nc_gene], ]), na.rm = T)
			numerator = sum(Model$y[g, ] * numerator, na.rm = T)
			numerator = mean_xi / (sd_xi ^ 2) + numerator
			# calculating denominator
			denominator = sum(colSums(as.matrix(tau[gihj, ]), na.rm = T) * Model$y2[g, ], na.rm = T)
			denominator = 1 / (sd_xi ^ 2) + denominator
			x_seq = numerator / denominator
			x2_seq = x_seq ^ 2 + 1 / denominator
			return(list(x_seq = x_seq, x2_seq = x2_seq))
		}
	
	xx_loop <-
		function(i,
				 Model,
				 LFC,
				 x_seq1_genes,
				 x_seq2_genes,
				 pattern_join,
				 tau,
				 mean_xx,
				 sd_xx) {
			# Getting individual mean_xxij and sd_xxij values for each guide pair ij
			mean_xxij = mean_xx[i]
			sd_xxij = sd_xx[i]
			
			# Get genes gh and guides ij
			gi = Model$hashes_x$paired_guide[i, 1]
			hj = Model$hashes_x$paired_guide[i, 2]
			g = x_seq1_genes[gi]
			h = x_seq2_genes[hj]
			gh = paste(sort(c(g, h)), collapse = pattern_join)
			
			numerator = Model$s[gh, ] * tau[i, ] * (LFC[i, ] - Model$x_seq1[gi] *
														Model$y[g, ] - Model$x_seq2[hj] * Model$y[h, ])
			numerator = mean_xxij / (sd_xxij ^ 2) + sum(numerator, na.rm = T)
			denominator = 1 / (sd_xxij ^ 2) + sum(Model$s2[gh, ] * tau[i, ], na.rm = T)
			xx = numerator / denominator
			xx2 = xx ^ 2 + 1 / denominator
			return(list(xx = xx, xx2 = xx2))
		}
	if (verbose) {
		message("\tUpdating x_seq1:")
	}
	x1_res <- pbmclapply(
		X = names(Model$x_seq1),
		FUN = x1_loop,
		Model = Model,
		guide2gene = guide2gene,
		LFC = LFC,
		tau = tau,
		x_seq1_genes = x_seq1_genes,
		pattern_join = pattern_join,
		nc_gene = nc_gene,
		mean_x = mean_x,
		sd_x = sd_x,
		mc.cores = cores,
		ignore.interactive = verbose
	)
	
	# Note: use [] to coerce output to vector
	Model$x_seq1[] <- lapply(x1_res, magrittr::extract, "x_seq") %>%
		unlist(recursive = T, use.names = F)
	Model$x2_seq1[] <- lapply(x1_res, magrittr::extract, "x2_seq") %>%
		unlist(recursive = T, use.names = F)
	
	if (verbose) {
		message("\tUpdating x_seq2:")
	}
	x2_res <- pbmclapply(
		X = names(Model$x_seq2),
		FUN = x2_loop,
		Model = Model,
		guide2gene = guide2gene,
		LFC = LFC,
		tau = tau,
		x_seq2_genes = x_seq2_genes,
		pattern_join = pattern_join,
		nc_gene = nc_gene,
		mean_x = mean_x,
		sd_x = sd_x,
		mc.cores = cores,
		ignore.interactive = verbose
	)
	
	Model$x_seq2[] <- lapply(x2_res, magrittr::extract, "x_seq") %>%
		unlist(recursive = T, use.names = F)
	Model$x2_seq2[] <- lapply(x2_res, magrittr::extract, "x2_seq") %>%
		unlist(recursive = T, use.names = F)
	
	if (verbose) {
		message("\tUpdating xx:")
	}
	xx_res <- pbmclapply(
		names(Model$xx),
		FUN = xx_loop,
		Model = Model,
		LFC = LFC,
		x_seq1_genes = x_seq1_genes,
		x_seq2_genes = x_seq2_genes,
		pattern_join = pattern_join,
		tau = tau,
		mean_xx = mean_xx,
		sd_xx = sd_xx,
		mc.cores = cores,
		ignore.interactive = verbose
	)
	
	Model$xx[] <- lapply(xx_res, extract, "xx") %>%
		unlist(use.names = F)
	Model$xx2[] <- lapply(xx_res, extract, "xx2") %>%
		unlist(use.names = F)
	
	# output
	if (verbose) {
		tend = Sys.time()
		tdiff = difftime(tend, tstart)
		message("\tCompleted update of x.")
		message("\tTime to completion: ",
				round(tdiff, digits = 3),
				' ',
				units(tdiff))
	}
	return(Model)
}
