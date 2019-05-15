#' initialize_tau
#' @description Initialize all tau values based on the observed replicate variance.
#'
#' @param Model an object of class gemini.model
#' @param CONSTANT a numeric indicating a constant value that shifts counts to reduce outliers (default = 32).
#' @param prior_shape shape parameter of Gamma distribution used to model the variation in the data in \code{Input}. If single numeric value, then shape parameters for all samples are assumed equal. Otherwise, a named numeric vector of shape parameters the same length as the number of samples (excluding early time point).
#' @param window numeric if window smoothing should be done on initialized tau values, otherwise NULL (default) for no window smoothing
#' @param monotonize logical specifying whether the variance should be monotonically increasing (default FALSE)
#' @param verbose default FALSE
#'
#' @return a Model object of class gemini.model including new slots for alpha and beta values
#' 
#' @importFrom stats sd
#' @export
#' @examples
#' data("Model", package = "gemini")
#' Model <- initialize_tau(Model, CONSTANT = 32, prior_shape = 0.5)
initialize_tau <- function(Model,
						   CONSTANT = 32,
						   prior_shape = 0.5,
						   window = NULL,
						   monotonize = FALSE,
						   verbose = FALSE) {
	# Check input
	stopifnot("gemini.model" %in% class(Model))
	
	# User message
	if (verbose)
		message("Initializing tau")
	
	Input <- Model$Input
	sample2replicate <- Input[[Model$replicate.map]]
	
	# cell line names
	sample_early = unique(sample2replicate[sample2replicate[, "TP"] == "ETP", Model$sample.column.name])
	sample_late = unique(sample2replicate[sample2replicate[, "TP"] == "LTP", Model$sample.column.name])
	
	# Pass to appropriate helper function, given the number of ETP samples
	if (length(intersect(sample_early, sample_late))==0) {
		Model <-
			initialize_tau_single_etp(
				Model,
				CONSTANT,
				prior_shape,
				window,
				monotonize,
				verbose
			)
	} else{
		Model <-
			initialize_tau_multi_etp(
				Model,
				CONSTANT,
				prior_shape,
				window,
				monotonize,
				verbose
			)
	}
	return(Model)
}

initialize_tau_multi_etp <- function(Model,
									 CONSTANT = 32,
									 prior_shape = 0.5,
									 window = NULL,
									 monotonize = FALSE,
									 verbose = FALSE) {
	Input <- Model$Input
	sample2replicate <- Input[[Model$replicate.map]]
	
	# cell line names
	sample_early = unique(sample2replicate[sample2replicate[, "TP"] == "ETP", Model$sample.column.name])
	sample_late = unique(sample2replicate[sample2replicate[, "TP"] == "LTP", Model$sample.column.name])
	
	# prior shape as vector or constant?
	if (length(prior_shape) == 1) {
		prior_shape <- rep(prior_shape, length(sample_late))
		names(prior_shape) <- sample_late
	} else if (length(prior_shape) != length(sample_late)) {
		stop(
			"Prior shape not defined for each sample! Please ensure that length of prior shape is the same as the number of LTP samples."
		)
	}
	if(!all(sample_late %in% names(prior_shape))){
		stop("Prior shape not defined for each LTP sample! Please ensure that length of prior shape is the same as the number of LTP samples and that names of prior_shape match sample names.")
	}
	Model$prior_shape <- prior_shape
	
	
	# define mean and sd matrix
	meanmat = sdmat = matrix(
		nrow = nrow(Input$counts),
		ncol = length(sample_late),
		dimnames = list(rownames(Input$counts), sample_late)
	)
	
	meanmat_etp = sdmat_etp = matrix(
		nrow = nrow(Input$counts),
		ncol = length(sample_early),
		dimnames = list(rownames(Input$counts), sample_early)
	)
	
	# average counts across all samples
	countavr = sum(Input$counts, na.rm = TRUE) / (ncol(Input$counts))
	
	# mean and sd calculation
	for (i in sample_early) {
		col.ind = sample2replicate[sample2replicate[, Model$sample.column.name] == i &
								   	sample2replicate[, "TP"] == "ETP", 1] # replicate names for cell line i
		count_replicates = as.matrix(Input$counts[, col.ind])
		count_replicates = apply(count_replicates, 2, function(x)
			x / sum(x, na.rm = TRUE)) * countavr # normalize samples by total counts
		count_replicates_mednorm = apply(count_replicates, 2, function(x)
			.median_normalize(x, CONSTANT)) # median-normalize samples
		meanmat_etp[, i] = rowMeans(count_replicates_mednorm, na.rm = TRUE) # calculate mean across rows
		sdmat_etp[, i] = apply(count_replicates_mednorm, 1, function(x)
			sd(x, na.rm = TRUE)) # calculate sd across rows
	}
	
	for (i in sample_late) {
		col.ind = sample2replicate[sample2replicate[, Model$sample.column.name] == i &
								   	sample2replicate[, "TP"] == "LTP", 1] # replicate names for cell line i
		count_replicates = as.matrix(Input$counts[, col.ind])
		count_replicates = apply(count_replicates, 2, function(x)
			x / sum(x, na.rm = TRUE)) * countavr # normalize samples by total counts
		count_replicates_mednorm = apply(count_replicates, 2, function(x)
			.median_normalize(x, CONSTANT)) # median-normalize samples
		meanmat[, i] = rowMeans(count_replicates_mednorm, na.rm = TRUE) # calculate mean
		sdmat[, i] = apply(count_replicates_mednorm, 1, function(x)
			sd(x, na.rm = TRUE)) # calculate sd
	}
	
	# window smoothing and monotonizing
	if (!is.null(window)) {
		for (i in seq_len(length(sample_early))) {
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
		for (i in seq_len(length(sample_late))) {
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
		sdmat[, colSums(!is.na(sdmat)) == 0] = median(cbind(sdmat, sdmat_etp), na.rm = TRUE)
		sdmat_etp[, colSums(!is.na(sdmat_etp)) == 0] = median(cbind(sdmat, sdmat_etp), na.rm = TRUE)
		
		# median uncertainty within one cell line if no replicate available for a guide
		sdmat = apply(sdmat, 2, function(x) {
			x[is.na(x)] = median(x, na.rm = TRUE)
			return(x)
		})
		sdmat_etp = apply(sdmat_etp, 2, function(x) {
			x[is.na(x)] = median(x, na.rm = TRUE)
			return(x)
		})
		
		beta <-
			matrix(
				nrow = nrow(sdmat),
				ncol = length(sample_late),
				dimnames = list(rownames(sdmat), sample_late)
			)
		for (i in sample_late) {
			beta[, i] = prior_shape[i] * sdmat[, i] ^ 2 + sdmat_etp[, sample_early == i] ^
				2 + 2 * 1e-2
		}
	} else {
		# if no replicate available for any cell line
		beta = matrix(
			0.2 * prior_shape[sample_late],
			nrow = nrow(sdmat),
			ncol = ncol(sdmat),
			dimnames = list(rownames(sdmat), sample_late),
			byrow = TRUE
		) # assume beta = 0.2*prior_shape, leading to the variance beta/alpha = 0.2 for LFCs
	}
	
	alpha = matrix(
		prior_shape[sample_late],
		nrow = nrow(beta),
		ncol = ncol(beta),
		dimnames = list(rownames(beta), colnames(beta)),
		byrow = TRUE
	)
	
	# output
	Model$alpha = alpha
	Model$beta = Model$beta_prior = beta
	return(Model)
}

initialize_tau_single_etp <- function(Model,
									  CONSTANT = 32,
									  prior_shape = 0.5,
									  window = NULL,
									  monotonize = FALSE,
									  verbose = FALSE) {
	Input <- Model$Input
	sample2replicate <- Input[[Model$replicate.map]]
	
	# cell line names
	sample_early = unique(sample2replicate[sample2replicate[, "TP"] == "ETP", Model$sample.column.name])
	sample = unique(sample2replicate[,Model$sample.column.name])
	
	# prior shape as vector or constant?
	if (length(prior_shape) == 1) {
		prior_shape <-
			rep(prior_shape, length(unique(sample2replicate$samplename[sample2replicate[, "TP"] == 'LTP'])))
		names(prior_shape) <- unique(sample2replicate$samplename[sample2replicate[, "TP"] == 'LTP'])
	} else if (length(prior_shape) != length(unique(sample2replicate[sample2replicate[, "TP"] ==
																	 "LTP", Model$sample.column.name]))) {
		stop(
			"Prior shape not defined for each LTP sample! Please ensure that length of prior shape is the same as the number of LTP samples, or that one prior shape is defined for all."
		)
	}
	if(!all(unique(sample2replicate[sample2replicate[, "TP"] == "LTP", Model$sample.column.name]) %in% names(prior_shape))){
		stop("Prior shape not defined for each LTP sample! Please ensure that length of prior shape is the same as the number of LTP samples and that names of prior_shape match sample names.")
	}
	Model$prior_shape <- prior_shape
	
	# define mean and sd matrix
	meanmat = sdmat = matrix(
		nrow = nrow(Input$counts),
		ncol = length(sample),
		dimnames = list(rownames(Input$counts), sample)
	)
	countavr = sum(Input$counts, na.rm = TRUE) / (ncol(Input$counts)) # average counts per replicate
	
	# mean and sd calculation
	for (i in sample) {
		col.ind = sample2replicate[sample2replicate[, Model$sample.column.name] == i, 1] # replicate names for sample i
		count_replicates = as.matrix(Input$counts[, col.ind])
		count_replicates = apply(count_replicates, 2, function(x)
			x / sum(x, na.rm = TRUE)) * countavr # normalize samples by total counts
		count_replicates_mednorm = apply(count_replicates, 2, function(x)
			.median_normalize(x, CONSTANT)) # median-normalize samples
		meanmat[, i] = rowMeans(count_replicates_mednorm, na.rm = TRUE) # calculate mean
		sdmat[, i] = apply(count_replicates_mednorm, 1, function(x)
			sd(x, na.rm = TRUE)) # calculate sd
	}
	
	# window smoothing and monotonizing
	if (!is.null(window)) {
		for (i in seq_len(length(sample))) {
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
		sdmat[, colSums(!is.na(sdmat)) == 0] = median(sdmat, na.rm = TRUE)
		
		# median uncertainty within one cell line if no replicate available for a guide
		sdmat = apply(sdmat, 2, function(x) {
			x[is.na(x)] = median(x, na.rm = TRUE)
			return(x)
		})
		
		beta = t(t(sdmat[,!(sample_early == colnames(sdmat))] ^ 2 + sdmat[, sample_early] ^
				   	2 + 2 * 1e-2) * prior_shape[sample[!(sample_early == colnames(sdmat))]]) %>%
			as.matrix()  %>%
			set_colnames(sample[!(sample_early == colnames(sdmat))]) %>% set_rownames(rownames(sdmat))
	} else {
		# if no replicate available for any cell line
		beta = matrix(
			0.2 * prior_shape[sample[!(sample_early == colnames(sdmat))]],
			nrow = dim(sdmat)[1],
			ncol = sum(!(sample_early == colnames(sdmat))),
			dimnames = list(rownames(sdmat), sample[!(sample_early == colnames(sdmat))]),
			byrow = TRUE
		) # assume beta = 0.2*prior_shape, leading to the variance beta/alpha = 0.2 for LFCs
	}
	
	alpha = matrix(
		prior_shape[sample[!(sample_early == colnames(sdmat))]],
		nrow = dim(beta)[1],
		ncol = dim(beta)[2],
		dimnames = list(rownames(beta), colnames(beta)),
		byrow = TRUE
	)
	
	# output
	Model$alpha = alpha
	Model$beta = Model$beta_prior = beta
	return(Model)
}
