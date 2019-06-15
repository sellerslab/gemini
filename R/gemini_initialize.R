#' gemini_initialize
#' @description Creates a gemini.model object given Input data and initialization parameters.
#'
#' @param Input an object of class gemini.input
#' @param guide.pair.annot the name of an object within \code{Input} containing guide to gene annotations (Default = "guide.pair.annot")
#' @param replicate.map the name of an object within \code{Input} containing replicate and sample mappings
#' (Default = "replicate.map")
#' @param sample.column.name a character or integer indicating which column of \code{Input$replicate.map} describes the samples.
#' @param LFC.name a character indicating an object within Input to treat as LFC. By default, "LFC" is used,
#' which is the output of \code{\link[gemini]{gemini_calculate_lfc}}.
#' @param nc_gene a character naming the gene to use as a negative control. See details for more.
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
#' @param cores numeric indicating the number of cores to use. See details in \code{\link[gemini]{gemini_parallelization}}.
#' @param save a character (file path) or logical indicating whether the initialized model should be saved.
#'
#' @details
#' guide.pair.annot is created with the following format:
#' \itemize{
#'  \item 1st column contains guide pairs, joined by \code{pattern_split}.
#'  \item 2nd column contains the name of gene mapping to the first guide in the 1st column.
#'  \item 3rd column contains the name of gene mapping to the second guide in the 1st column.
#' }
#' Additional columns may be appended to guide.pair.annot, but the first three columns must be defined as
#' above.
#'
#' \strong{Use of negative control}{
#' As described in Zamanighomi et al., it is highly recommended to specify a negative control. In the event that no negative control is specified, GEMINI will use the median LFC of all guide pairs targeting a gene to initialize that gene's effect. While this may be reasonable in large all-by-all screens, it is \strong{not recommended} in smaller screens or some-by-some screens. As a result, when possible, be sure to specify a negative control.
#' }
#'
#' @return a Model object of class gemini.model
#' @export
#' 
#' @examples 
#' data("Input", package = "gemini")
#' Model <- gemini_initialize(Input, nc_gene = "CD81")
#' 
gemini_initialize <- function(Input,
                              guide.pair.annot = "guide.pair.annot",
                              replicate.map = "replicate.map",
                              sample.column.name = "samplename",
                              LFC.name = "LFC",
                              nc_gene,
                              CONSTANT = 32,
                              concordance = 1,
                              prior_shape = 0.5,
                              pattern_join = ";",
                              pattern_split = ";",
                              window = NULL,
                              monotonize = FALSE,
                              verbose = FALSE,
                              cores = 1,
                              save = NULL) {
    # Check Input
    stopifnot("gemini.input" %in% class(Input))
    
    # Check parallelization - always leave 1 core free
    if (cores >= parallel::detectCores()) {
        cores = parallel::detectCores() - 1
    }
    
    if (!"LFC" %in% names(Input)) {
        stop("No LFC found! Did you run gemini_calculate_lfc?")
    }
    
    # Create an empty model object
    Model <- list()
    class(Model) <- c(class(Model), "gemini.model")
    
    # Initialize all Input values
    Model$Input <- Input
    Model$guide.pair.annot = guide.pair.annot
    Model$replicate.map = replicate.map
    Model$sample.column.name = sample.column.name
    Model$nc_gene = nc_gene
    Model$LFC.name = LFC.name
    Model$pattern_join = pattern_join
    Model$pattern_split = pattern_split
    
    # Run initialization functions
    Model %<>%
        initialize_x(concordance = concordance,
                     cores = cores,
                     verbose = verbose) %>%
        initialize_y(verbose = verbose,
                     cores = cores) %>%
        initialize_s(verbose = verbose,
                     cores = cores) %>%
        initialize_tau(
            CONSTANT = CONSTANT,
            prior_shape = prior_shape,
            window = window,
            monotonize = monotonize,
            verbose = verbose
        )
    
    # initialize mae
    Model$mae <- numeric()
    
    # mae calculated using initial values
    Model <- update_mae(Model = Model,
                        verbose = verbose)
    
    if (verbose)
        message("Model initialized.")
    
    if (!is.null(save)) {
        # Save initialized model
        if (!is.character(save) & length(save) == 1) {
            save = "gemini_initialized_model"
        }
        saveRDS(Model, file = save)
    }
    
    return(Model)
}
