#' initialize_tau
#' @description Initialize all tau values based on the observed replicate variance.
#'
#' @param Model an object of class gemini.model
#' @param Input an object of class gemini.input
#' @param replicate.map an object within \code{Input} containing replicate and sample mappings
#' @param sample.column.name a character or integer indicating which column of \code{Input$replicate.map} describes the samples.
#' @param CONSTANT a numeric indicating a constant value that shifts counts to reduce outliers (default = 32).
#' @param prior_shape shape parameter of Gamma distribution used to model the variation in the data in \code{Input}
#' @param window numeric if window smoothing should be done on initialized tau values, otherwise NULL (default) for no window smoothing
#' @param monotonize logical specifying whether the variance should be monotonically increasing (default FALSE)
#' @param verbose default FALSE
#'
#' @return a Model object of class gemini.model including new slots for alpha and beta values
#' @export
initialize_tau <- function(Model,
						   Input,
						   replicate.map = "replicate.map",
						   sample.column.name = "samplename",
						   CONSTANT = 32,
						   prior_shape = 0.5,
						   window = NULL,
						   monotonize = FALSE, verbose = F){
	# Check input
	stopifnot("gemini.model" %in% class(Model))
	
	# User message
	if(verbose) message("Initializing tau")

	sample2replicate <- Input[[replicate.map]]
	
	# cell line names
	sample_early = unique(sample2replicate[sample2replicate[,"TP"]=="ETP", sample.column.name])
	cell_lines = unique(sample2replicate[sample2replicate[,"TP"]=="LTP", sample.column.name])

	# Pass to appropriate helper function, given the number of ETP samples
	if(length(sample_early)==1){
		Model <- initialize_tau_single_etp(Model, Input, replicate.map, sample.column.name, CONSTANT, prior_shape, window, monotonize, verbose)
	}else{
		Model <- initialize_tau_multi_etp(Model, Input, replicate.map, sample.column.name, CONSTANT, prior_shape, window, monotonize, verbose)
	}
	return(Model)
}

initialize_tau_multi_etp <- function(Model,
									 Input,
									 replicate.map = "replicate.map",
									 sample.column.name = "samplename",
									 CONSTANT = 32,
									 prior_shape = 0.5,
									 window = NULL,
									 monotonize = FALSE, verbose = F){
	sample2replicate <- Input[[replicate.map]]
	
	# cell line names
	sample_early = unique(sample2replicate[sample2replicate[,"TP"]=="ETP", sample.column.name])
	cell_lines = unique(sample2replicate[sample2replicate[,"TP"]=="LTP", sample.column.name])
	
	# prior shape as vector or constant?
	if(length(prior_shape)==1){
		prior_shape <- rep(prior_shape, length(cell_lines))
	}else if(length(prior_shape)!=length(cell_lines)){
		stop("Prior shape not defined for each LTP cell line! Please ensure that length of prior shape is the same as the number of LTP samples.")
	}
	Model$prior_shape <- prior_shape
	
	
	# define mean and sd matrix
	meanmat = sdmat = matrix(
		nrow = nrow(Input$counts),
		ncol = length(cell_lines),
		dimnames = list(rownames(Input$counts), cell_lines)
	)
	
	meanmat_etp = sdmat_etp = matrix(
		nrow = nrow(Input$counts),
		ncol = length(sample_early),
		dimnames = list(rownames(Input$counts), sample_early)
	)
	
	# average counts across all samples
	countavr = sum(Input$counts, na.rm = T) / (ncol(Input$counts)) 
	
	# mean and sd calculation
	for (i in sample_early) {
		col.ind = sample2replicate[sample2replicate[, sample.column.name] == i & sample2replicate[, "TP"] == "ETP", 1] # replicate names for cell line i
		count_replicates = as.matrix(Input$counts[, col.ind])
		count_replicates = apply(count_replicates, 2, function(x)
			x / sum(x, na.rm = T)) * countavr # normalize samples by total counts
		count_replicates_mednorm = apply(count_replicates, 2, function(x)
			.median_normalize(x, CONSTANT)) # median-normalize samples
		meanmat_etp[, i] = rowMeans(count_replicates_mednorm, na.rm = T) # calculate mean across rows
		sdmat_etp[, i] = apply(count_replicates_mednorm, 1, function(x)
			sd(x, na.rm = T)) # calculate sd across rows
	}
	
	for (i in cell_lines) {
		col.ind = sample2replicate[sample2replicate[, sample.column.name] == i & sample2replicate[, "TP"] == "LTP", 1] # replicate names for cell line i
		count_replicates = as.matrix(Input$counts[, col.ind])
		count_replicates = apply(count_replicates, 2, function(x)
			x / sum(x, na.rm = T)) * countavr # normalize samples by total counts
		count_replicates_mednorm = apply(count_replicates, 2, function(x)
			.median_normalize(x, CONSTANT)) # median-normalize samples
		meanmat[, i] = rowMeans(count_replicates_mednorm, na.rm = T) # calculate mean
		sdmat[, i] = apply(count_replicates_mednorm, 1, function(x)
			sd(x, na.rm = T)) # calculate sd
	}
	
	# window smoothing and monotonizing
	if (!is.null(window)) {
		for (i in 1:length(sample_early)) {
			I = order(meanmat_etp[, i], decreasing = FALSE)
			vars = sdmat_etp[I, i] ^ 2
			# window smoothing
			vars = .window_smooth(vars, window)
			# monotonizing
			if (monotonize) {
				vars = .monotonize(vars)
			}
			sdmat_etp[I, i] = vars ^ 0.5
		}
		for (i in 1:length(cell_lines)) {
			I = order(meanmat[, i], decreasing = FALSE)
			vars = sdmat[I, i] ^ 2
			# window smoothing
			vars = .window_smooth(vars, window)
			# monotonizing
			if (monotonize) {
				vars = .monotonize(vars)
			}
			sdmat[I, i] = vars ^ 0.5
		}
	}
	
	# sd calculation when one replicate available
	if (any(!is.na(sdmat)) | any(!is.na(sdmat_etp))) {
		# at least two replicates available for minimum one cell line
		# median uncertainty across all cell lines if no replicate available for a cell line
		sdmat[, colSums(!is.na(sdmat)) == 0] = median(cbind(sdmat, sdmat_etp), na.rm = T)
		sdmat_etp[, colSums(!is.na(sdmat_etp)) == 0] = median(cbind(sdmat, sdmat_etp), na.rm = T)
		
		# median uncertainty within one cell line if no replicate available for a guide
		sdmat = apply(sdmat, 2, function(x) {
			x[is.na(x)] = median(x, na.rm = T)
			return(x)
		})
		sdmat_etp = apply(sdmat_etp, 2, function(x) {
			x[is.na(x)] = median(x, na.rm = T)
			return(x)
		})
		
		beta <- matrix(nrow = nrow(sdmat), ncol = length(cell_lines), dimnames = list(rownames(sdmat), cell_lines))
		for(i in cell_lines){
			beta[,i] = prior_shape[i] * sdmat[,i]^2 + sdmat_etp[,sample_early == i]^2 + 2 * 1e-2
		}
	} else {
		# if no replicate available for any cell line
		beta = matrix(
			0.2 * prior_shape,
			nrow = nrow(sdmat),
			ncol = ncol(sdmat),
			dimnames = list(rownames(sdmat), cell_lines),
			byrow = T
		) # assume beta = 0.2*prior_shape, leading to the variance beta/alpha = 0.2 for LFCs
	}
	
	alpha = matrix(
		prior_shape,
		nrow = nrow(beta),
		ncol = ncol(beta),
		dimnames = list(rownames(beta), colnames(beta)),
		byrow = T
	)
	
	# output
	Model$alpha = alpha
	Model$beta = Model$beta_prior = beta
	return(Model)
}

initialize_tau_single_etp <- function(Model,
									  Input,
									  replicate.map = "replicate.map",
									  sample.column.name = "samplename",
									  CONSTANT = 32,
									  prior_shape = 0.5,
									  window = NULL,
									  monotonize = FALSE, verbose = F){
	sample2replicate <- Input[[replicate.map]]
	
	# cell line names
	sample_early = unique(sample2replicate[sample2replicate[,"TP"]=="ETP", sample.column.name])
	cell_lines = unique(sample2replicate[, sample.column.name])
	
	# prior shape as vector or constant?
	if(length(prior_shape)==1){
		prior_shape <- rep(prior_shape, length(unique(sample2replicate$samplename[sample2replicate[,"TP"]=='LTP'])))
	}else if(length(prior_shape)!=length(unique(sample2replicate[sample2replicate[,"TP"]=="LTP", sample.column.name]))){
		stop("Prior shape not defined for each LTP cell line! Please ensure that length of prior shape is the same as the number of LTP samples.")
	}
	Model$prior_shape <- prior_shape
	
	# define mean and sd matrix
	meanmat = sdmat = matrix(
		nrow = nrow(Input$counts),
		ncol = length(cell_lines),
		dimnames = list(rownames(Input$counts), cell_lines)
	)
	countavr = sum(Input$counts, na.rm = T) / (ncol(Input$counts)) # average counts across all samples
	
	# mean and sd calculation
	for (i in cell_lines) {
		col.ind = sample2replicate[sample2replicate[, sample.column.name] == i, 1] # replicate names for cell line i
		count_replicates = as.matrix(Input$counts[, col.ind])
		count_replicates = apply(count_replicates, 2, function(x)
			x / sum(x, na.rm = T)) * countavr # normalize samples by total counts
		count_replicates_mednorm = apply(count_replicates, 2, function(x)
			.median_normalize(x, CONSTANT)) # median-normalize samples
		meanmat[, i] = rowMeans(count_replicates_mednorm, na.rm = T) # calculate mean
		sdmat[, i] = apply(count_replicates_mednorm, 1, function(x)
			sd(x, na.rm = T)) # calculate sd
	}
	
	# window smoothing and monotonizing
	if (!is.null(window)) {
		for (i in 1:length(cell_lines)) {
			I = order(meanmat[, i], decreasing = FALSE)
			vars = sdmat[I, i] ^ 2
			# window smoothing
			vars = .window_smooth(vars, window)
			# monotonizing
			if (monotonize) {
				vars = .monotonize(vars)
			}
			sdmat[I, i] = vars ^ 0.5
		}
	}
	
	# sd calculation when one replicate available
	if (any(!is.na(sdmat))) {
		# at least two replicates available for minimum one cell line
		# median uncertainty across all cell lines if no replicate available for a cell line
		sdmat[, colSums(!is.na(sdmat)) == 0] = median(sdmat, na.rm = T)
		# median uncertainty within one cell line if no replicate available for a guide
		sdmat = apply(sdmat, 2, function(x) {
			x[is.na(x)] = median(x, na.rm = T)
			return(x)
		})
		
		beta = t(t(sdmat[, !(sample_early == colnames(sdmat))] ^ 2 + sdmat[, sample_early] ^
							  	2 + 2 * 1e-2) * prior_shape) %>%
			as.matrix()  %>%
			set_colnames(cell_lines[!(sample_early == colnames(sdmat))]) %>% set_rownames(rownames(sdmat))
	} else {
		# if no replicate available for any cell line
		beta = matrix(
			0.2 * prior_shape,
			nrow = dim(sdmat)[1],
			ncol = sum(!(sample_early == colnames(sdmat))),
			dimnames = list(rownames(sdmat), cell_lines[!(sample_early == colnames(sdmat))]),
			byrow = T
		) # assume beta = 0.2*prior_shape, leading to the variance beta/alpha = 0.2 for LFCs
	}
	
	alpha = matrix(
		prior_shape,
		nrow = dim(beta)[1],
		ncol = dim(beta)[2],
		dimnames = list(rownames(beta), colnames(beta)),
		byrow = T
	)
	
	# output
	Model$alpha = alpha
	Model$beta = Model$beta_prior = beta
	return(Model)
}