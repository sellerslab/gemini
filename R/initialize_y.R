#' initialize_y
#' @description Initialize all y values from guide pairs including a negative control.
#'
#' @param Model an object of class gemini.model
#' @param verbose default FALSE
#' @param cores a numeric indicating the number of cores to use.  See details in \code{\link[gemini]{gemini_parallelization}}. (default=1)
#'
#' @return a Model object of class gemini.model including new slots for y values
#' 
#' @examples 
#' data("Model", package = "gemini")
#' Model %<>% initialize_y()
#'
#' @export
initialize_y <- function(Model,
                         verbose = FALSE,
                         cores = 1) {
    # Check input
    stopifnot("gemini.model" %in% class(Model))
    if (verbose)
        message("Initializing y")
    
    # Create inputs
    Input <- Model$Input
    guide2gene <- Input[[Model$guide.pair.annot]]
    LFC <- Input[[Model$LFC.name]]
    
    # gene pairs containing nc_gene
    guide2nc_gene = guide2gene[guide2gene[, 2] %in% Model$nc_gene |
                                   guide2gene[, 3] %in% Model$nc_gene, ]
    
    # genes mapped to guides
    hash_nc = Sgene2Pguides_hash(guide2nc_gene, cores)
    hash_others = Sgene2Pguides_hash(guide2gene, cores)
    
    # gene list
    genes = unique(c(guide2gene[, 2], guide2gene[, 3]))
    
    # removing nc_gene from gene list
    genes = genes[!is.element(genes, Model$nc_gene)]
    
    # define matrix y
    y = matrix(
        nrow = length(genes),
        ncol = ncol(LFC),
        dimnames = list(genes, colnames(LFC))
    )
    
    # estimate y
    for (i in genes) {
        # all genes should be paired with nc_gene
        ind = hash_nc[[i]]
        if (length(ind) > 1) {
            y[i, ] = apply(as.matrix(LFC[ind, ]), 2, function(x)
                median(x, na.rm = TRUE))
        } else if (length(ind) == 1){
            y[i, ] = apply(as.matrix(t(LFC[ind, ])), 2, function(x)
                median(x, na.rm = TRUE))
        } else{
            # ... but if not, use the median LFC of all interactions involving the gene
            # NOTE: only a good assumption if you have an all-by-all screen!
            ind = hash_others[[i]]
            if (length(ind) > 1) {
                y[i, ] = apply(as.matrix(LFC[ind, ]), 2, function(x)
                    median(x, na.rm = TRUE))
            } else{
                y[i, ] = apply(as.matrix(t(LFC[ind, ])), 2, function(x)
                    median(x, na.rm = TRUE))
            }
        }
    }
    
    # output
    Model$y <- y
    Model$y2 <- y ^ 2
    
    return(Model)
}
