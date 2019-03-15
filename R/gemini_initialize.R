#' gemini_initialize
#' @description Creates a gemini.model object given Input data and initialization parameters.
#' 
#' @param Input an object of class gemini.input
#' @param guide.pair.annot the name of an object within \code{Input} containing guide to gene annotations (Default = "guide.pair.annot")
#' @param replicate.map the name of an object within \code{Input} containing replicate and sample mappings
#' (Default = "replicate.map")
#' @param sample.column.name a character or integer indicating which column of \code{Input$replicate.map} describes the samples.
#' @param nc_gene a character naming the gene to use as a negative control - REQUIRED
#' @param CONSTANT a numeric indicating a constant value that shifts counts to reduce outliers, see Details (default = 32).
#' @param concordance a numeric value to initialize x (default = 1)
#' @param prior_shape shape parameter of Gamma distribution used to model the variation in the 
#' data in \code{Input} (default = 0.5)
#' @param pattern_join a character to join the gene combinations found in 
#' \code{guide.pair.annot}. Default ';'
#' @param pattern_split a character to split the guide combinations found in 
#' \code{guide.pair.annot}. Default ';'
#' @param window a numeric if window smoothing should be done on initialized tau values, 
#' otherwise NULL (default) for no window smoothing
#' @param monotonize a logical specifying whether the variance should be monotonically increasing. (default FALSE)
#' @param verbose default FALSE
#' @param cores numeric indicating the number of cores to use. See details in \code{\link[gemini]{gemini_parallel}}.
#' @param save a character (file path) or logical indicating whether the initialized model should be saved.
#'
#' @return a Model object of class gemini.model
#'
#' @export
gemini_initialize <- function(Input,
							  guide.pair.annot = "guide.pair.annot",
							  replicate.map = "replicate.map",
							  sample.column.name = "samplename",
							  nc_gene, CONSTANT = 32, concordance = 1, prior_shape = 0.5,
							  pattern_join = ";", pattern_split = ";",
							  window = NULL, monotonize = FALSE, verbose = F,
							  cores = 1, save = NULL){
	# Check Input
	stopifnot("gemini.input" %in% class(Input))
  
  if(!"LFC" %in% names(Input)){
    stop("No LFC found! Did you run gemini_calculate_lfc?")
  }
	
	# Create an empty model object
	Model <- list()
	class(Model) <- c(class(Model), "gemini.model")
	
	# Initialize all Input values
	LFC = "LFC"
	Model$Input <- Input
	Model$guide.pair.annot = guide.pair.annot
	Model$replicate.map = replicate.map
	Model$nc_gene = nc_gene
	Model$LFC = LFC
	Model$pattern_join = pattern_join
	Model$pattern_split = pattern_split

	# Run initialization functions
	Model %<>%
		initialize_x(Input = Input,
					 LFC = LFC,
					 guide.pair.annot = guide.pair.annot,
					 concordance = concordance,
					 pattern_split = pattern_split, cores = cores,
					 nc_gene = nc_gene, verbose = verbose) %>%
		initialize_y(Input = Input, LFC = LFC,
					 guide.pair.annot = guide.pair.annot, 
					 nc_gene = nc_gene, verbose = verbose, 
					 cores = cores) %>%
		initialize_s(Input = Input,
					 LFC = LFC, guide.pair.annot = guide.pair.annot,
					 nc_gene = nc_gene, pattern_join = pattern_join, 
					 verbose = verbose, cores = cores) %>%
		initialize_tau(Input = Input,
					   replicate.map = replicate.map,
					   sample.column.name = sample.column.name,
					   CONSTANT = CONSTANT, prior_shape = prior_shape,
					   window = window, monotonize = monotonize, 
					   verbose = verbose)
	
	# initialize mae
	Model$mae <- numeric()
	
	# mae calculated using initial values
	Model <- update_mae(Model = Model, Input = Input, LFC = LFC, nc_gene = nc_gene, guide.pair.annot = guide.pair.annot, pattern_join = pattern_join, verbose = verbose)
	
	if(verbose) message("Model initialized.")
	
	if(!is.null(save)){
		# Save initialized model
		if(!is.character(save) & length(save)==1){
			save = "gemini_initialized_model"
		}
		saveRDS(Model, file = save)
	}
	
	return(Model)
}
