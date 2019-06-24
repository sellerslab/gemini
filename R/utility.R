#' Compound pipe
#' 
#' @name %<>%
#' @importFrom magrittr %<>%
#' @rdname compound
#' @param lhs,rhs An object modified in place, and a function to apply to it
#' @export
#' @examples 
#' a <- 1:5
#' a %<>% sum() # a is modified in place
#' 
NULL

#' Pipe
#' 
#' @importFrom magrittr %>%
#' @name %>%
#' @rdname pipe
#' @export
#' @param lhs,rhs An object, and a function to apply to it
#' @examples
#' a <- 1:5
#' b <- a %>% sum()
NULL

#' Median normalization
#' @name .median_normalize
#' @param counts a counts matrix derived from Input
#' @param CONSTANT a pseudo-count to use (default = 32)
#' @importFrom stats median
#' @return median-normalized counts object
#' @examples
#' \dontrun{
#' .median_normalize(rnorm(n = 1000))
#' }
.median_normalize <- function(counts, CONSTANT = 32) {
    m = log2(counts + CONSTANT) - median(log2(counts + CONSTANT), na.rm = TRUE)
    
    # output
    return(m)
}

.window_smooth <- function(x, window = 800) {
    #flog.debug("Smoothing")
    N = length(x)
    m = rep(0, N) # mean variance in window
    m[seq_len(window + 1)] = mean(x[seq_len(2 * window + 1)]) # left shoulder; assume constant start to avoid high variance
    m[(N - window):N] = mean(x[(N - 2 * window):N]) # right shoulder; assume constant start to avoid high variance
    for (k in (window + 2):(N - window - 1)) {
        m[k] = m[k - 1] + (x[k + window + 1] - x[k - window]) / (2 * window + 1)
    } # linear time average
    #for (k in (window+2):(N-window-1)) { m[k] = mean(x[(k-window):(k+window+1)])} # linear time average
    return(m)
}

# monotonizing
.monotonize <- function(x) {
    #flog.debug("Monotonizing")
    N = length(x)
    for (i in seq_len(N)) {
        if (!is.nan(x[N - i + 1])) {
            x[N - i] = max(x[N - i + 1], x[N - i])
        }
    }
    return(x)
}

### Hash functions:
#' Gene hashing
#' @name Sgene2Pguides_hash
#' @param guide2gene derived from Input object
#' @param cores number of cores to use (default 1)
#' @importFrom parallel mclapply
#' @return Gene hash/list
#' @examples 
#' \dontrun{
#' #' data("Input", package = "gemini")
#' Sgene2Pguides_hash(Input$guide.pair.annot[1:10,])
#' }

Sgene2Pguides_hash <- function(guide2gene, cores = 1) {
    # Single gene to pair guide hash:
    # First, take all unique genes involved in screen
    genes = unique(c(guide2gene[, 2], guide2gene[, 3]))
    
    # Now identify paired guides associated with each individual gene
    hash = parallel::mclapply(genes, function(x) {
        guide2gene[(guide2gene[, 2] == x | guide2gene[, 3] == x), 1]
    }, mc.cores = cores) %>% magrittr::set_names(genes)
    
    # output
    return(hash)
}

#' Guide hashing
#' @name Sguide2Pguides_hash
#' @param guide2gene derived from Input object
#' @param split character to split guides
#' @param cores number of cores to use (default 1)
#' 
#' @importFrom parallel mclapply
#' @importFrom magrittr set_rownames
#' @importFrom dplyr bind_rows
#' @return Guide hash/list
#' @examples
#' \dontrun{
#' data("Input", package = "gemini") 
#' Sguide2Pguides_hash(gemini::Input$guide.pair.annot[1:10,])
#' }

Sguide2Pguides_hash <- function(guide2gene,
                                split = ";",
                                cores = 1) {
    # split guide sequences according to the pattern "split"
    paired_guide = parallel::mclapply(guide2gene[, 1], function(x) {
        s = strsplit(x, split = split, fixed = TRUE)[[1]]
        data.frame(
            guide1.sequence = s[1],
            guide2.sequence = s[2],
            stringsAsFactors = FALSE
        )
    }, mc.cores = cores, mc.cleanup = TRUE) %>%
        dplyr::bind_rows() %>%
        magrittr::set_rownames(guide2gene[, 1])
    
    if (all(is.na(paired_guide$guide1.sequence)) |
        all(is.na(paired_guide$guide2.sequence))) {
        stop("NAs detected for guide splitting - is pattern_split correct?")
    }
    
    # guide.sequence mapped to the paired guides containing the single guide
    hash <- lapply(unique(c(paired_guide$guide1.sequence, paired_guide$guide2.sequence)), function(x){
        rownames(paired_guide)[paired_guide$guide1.sequence == x | paired_guide$guide2.sequence == x]
    }) %>% magrittr::set_names(unique(c(paired_guide$guide1.sequence, paired_guide$guide2.sequence)))

    return(list(
        hash = hash,
        paired_guide = paired_guide
    ))
}

