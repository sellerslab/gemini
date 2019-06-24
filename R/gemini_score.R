#' Scoring Combination Effect
#'
#' @description Score genetic interactions from a gemini.model and produce a gemini.score object, and generate p-values and FDRs if negative control
#' pairs (\code{nc_pairs}) are available.
#'
#' @param Model an object of class gemini.model
#' @param pc_gene a character vector of any length naming genes to use as positive controls
#' @param pc_weight a weighting applied to the positive control
#' (\code{pc_gene}) to filter genes whose individual phenotype is
#' more lethal than \eqn{pc_weight*y(pc_gene)}.
#' @param pc_threshold a numeric value to indicate the LFC corresponding to a positive control, if no pc_gene is specified.
#' @param nc_pairs a set of non-interacting gene pairs to define
#' statistical significance.
#'
#' @return An object of class gemini.score containing score values for
#' strong interactions and sensitive lethality and recovery, and if \code{nc_pairs} is specified, statistical significance for each scoring type.
#'
#' @importFrom magrittr set_rownames
#' @importFrom magrittr set_names
#' @importFrom magrittr subtract
#' @importFrom magrittr add
#' @importFrom mixtools normalmixEM
#' @importFrom dplyr bind_rows
#' @importFrom stats p.adjust
#' @importFrom stats pnorm
#' @importFrom stats quantile
#' @importFrom utils capture.output
#'
#' @export
#' 
#' @examples
#' data("Model", package = "gemini")
#' Score <- gemini_score(Model, pc_gene = "EEF2")
gemini_score <- function(Model,
                         pc_gene = NA,
                         pc_threshold = NULL,
                         pc_weight = 0.5,
                         nc_pairs = NA) {
    # Check input
    stopifnot("gemini.model" %in% class(Model))
    
    Score <- list()
    class(Score) <- c(class(Score), "gemini.score")
    
    ###### Paired genes mapped to their corresponding Single genes
    Pgenes2Sgenes <- lapply(rownames(Model$s), function(x) {
        s <- strsplit(x, split = Model$pattern_join, fixed = TRUE)[[1]]
        data.frame(
            gene1 = s[1],
            gene2 = s[2],
            stringsAsFactors = FALSE
        )
    }) %>%
        dplyr::bind_rows() %>%
        magrittr::set_rownames(rownames(Model$s))
    
    ###### Score combination effects and calculate pvalue and fdr using nc_pairs
    Score$strong <- matrix(
        NA,
        nrow = nrow(Model$s),
        ncol = ncol(Model$s),
        dimnames = list(rownames(Model$s), colnames(Model$s))
    )
    Score$sensitive_lethality <- matrix(
        NA,
        nrow = nrow(Model$s),
        ncol = ncol(Model$s),
        dimnames = list(rownames(Model$s), colnames(Model$s))
    )
    Score$sensitive_recovery <- matrix(
        NA,
        nrow = nrow(Model$s),
        ncol = ncol(Model$s),
        dimnames = list(rownames(Model$s), colnames(Model$s))
    )
    Score$pvalue_strong <- matrix(
        NA,
        nrow = nrow(Model$s),
        ncol = ncol(Model$s),
        dimnames = list(rownames(Model$s), colnames(Model$s))
    )
    Score$pvalue_sensitive_lethality <- matrix(
        NA,
        nrow = nrow(Model$s),
        ncol = ncol(Model$s),
        dimnames = list(rownames(Model$s), colnames(Model$s))
    )
    Score$pvalue_sensitive_recovery <- matrix(
        NA,
        nrow = nrow(Model$s),
        ncol = ncol(Model$s),
        dimnames = list(rownames(Model$s), colnames(Model$s))
    )
    Score$fdr_strong <- matrix(
        NA,
        nrow = nrow(Model$s),
        ncol = ncol(Model$s),
        dimnames = list(rownames(Model$s), colnames(Model$s))
    )
    Score$fdr_sensitive_lethality <- matrix(
        NA,
        nrow = nrow(Model$s),
        ncol = ncol(Model$s),
        dimnames = list(rownames(Model$s), colnames(Model$s))
    )
    Score$fdr_sensitive_recovery <- matrix(
        NA,
        nrow = nrow(Model$s),
        ncol = ncol(Model$s),
        dimnames = list(rownames(Model$s), colnames(Model$s))
    )
    
    for (i in colnames(Model$s)) {
        # mean
        sgh_mean <- Model$s[, i]
        yg_mean <- Model$y[Pgenes2Sgenes[, 1], i]
        yh_mean <- Model$y[Pgenes2Sgenes[, 2], i]
        
        # sd
        sgh_sd <- (Model$s2[, i] - Model$s[, i] ^ 2) ^ 0.5
        yg_sd <- (Model$y2[Pgenes2Sgenes[, 1], i] - Model$y[Pgenes2Sgenes[, 1], i] ^
                     2) ^ 0.5
        yh_sd <- (Model$y2[Pgenes2Sgenes[, 2], i] - Model$y[Pgenes2Sgenes[, 2], i] ^
                     2) ^ 0.5
        
        # strong: Score synergies accroding to Score = abs(sgh_mean) - max(abs(yg_mean),abs(yh_mean)
        Score$strong[, i] <- abs(sgh_mean) - apply(cbind(abs(yg_mean), abs(yh_mean)), 1 , function(x)
            max(x, na.rm = TRUE)) %>%
            unlist() %>%
            magrittr::set_names(rownames(Model$s))
        # lethality: Score synergies accroding to Score = min(yg_mean,yh_mean) - yg_mean - yh_mean - sgh_mean
        Score$sensitive_lethality[, i] <- apply(cbind(yg_mean, yh_mean), 1 , function(x)
            min(x, na.rm = TRUE)) %>%
            unlist() %>%
            magrittr::set_names(rownames(Model$s)) %>%
            magrittr::subtract(sgh_mean) %>%
            magrittr::subtract(yg_mean) %>%
            magrittr::subtract(yh_mean)
        # recovery: Score synergies accroding to Score = yg_mean + yh_mean + sgh_mean - min(yg_mean,yh_mean)
        Score$sensitive_recovery[, i] <- sgh_mean %>%
            magrittr::add(yg_mean) %>%
            magrittr::add(yh_mean) %>%
            magrittr::subtract(apply(cbind(yg_mean, yh_mean), 1, function(x)
                min(x, na.rm = TRUE)) %>%
                    unlist() %>%
                    magrittr::set_names(rownames(Model$s)))
        
        # keep pairs if yg & yh > 0.5*ypc_gene
        if (sum(!is.na(pc_gene)) > 0 & is.null(pc_threshold)) {
            pc_gene <- pc_gene[!is.na(pc_gene)]
            threshold <- median(Model$y[pc_gene, i], na.rm = TRUE) * pc_weight
        } else if (sum(!is.na(pc_gene)) == 0 & !is.null(pc_threshold)) {
            threshold <- pc_threshold
        } else if (sum(!is.na(pc_gene)) > 0 & !is.null(pc_threshold)) {
            warning("Both pc_gene and pc_threshold specified - defaulting to pc_gene.")
            pc_gene <- pc_gene[!is.na(pc_gene)]
            threshold <- median(Model$y[pc_gene, i], na.rm = TRUE) * pc_weight
        }
        else {
            # If no threshold specified, use bottom 1-percentile of data
            threshold = as.numeric(quantile(Model$y[, i], probs = 0.01)) * pc_weight
        }
        remove_lethal = yg_mean < threshold | yh_mean < threshold
        remove_recovery = yg_mean > threshold & yh_mean > threshold
        Score$sensitive_lethality[remove_lethal , i] <- NA
        Score$sensitive_recovery[remove_recovery , i] <- NA
    } # end for
    
    ###### calculate pvalue and fdr using difference_constant or difference_quantile
    if (any(!is.na(nc_pairs))) {
        if (!requireNamespace("mixtools")) {
            stop("Please install the mixtools package: install.packages('mixtools').")
        }
        # construct negetaive model using mixture of two normal distributions
        # strong
        nc_pairs <- nc_pairs[!is.na(nc_pairs)]
        nc_values <- as.vector(Score$strong[nc_pairs, ])
        nc_values <- nc_values[!is.na(nc_values)]
        #cat("running EM algorithm on strong scoring \n")
        invisible(capture.output({
            nc_model_strong <- mixtools::normalmixEM(nc_values)
        }))
        # lethality
        nc_values <- as.vector(Score$sensitive_lethality[nc_pairs, ])
        nc_values <- nc_values[!is.na(nc_values)]
        #cat("running EM algorithm on lethality scoring \n")
        invisible(capture.output({
            nc_model_lethality <- mixtools::normalmixEM(nc_values)
        }))
        # recovery
        nc_values <- as.vector(Score$sensitive_recovery[nc_pairs, ])
        nc_values <- nc_values[!is.na(nc_values)]
        #cat("running EM algorithm on recovery ranking \n")
        invisible(capture.output({
            nc_model_recovery <- mixtools::normalmixEM(nc_values)
        }))
        
        # calculate p-values
        # strong
        Score$pvalue_strong <- nc_model_strong$lambda[1] * pnorm(
            Score$strong,
            mean = nc_model_strong$mu[1],
            sd = nc_model_strong$sigma[1],
            lower.tail = FALSE
        ) +
            nc_model_strong$lambda[2] * pnorm(
                Score$strong,
                mean = nc_model_strong$mu[2],
                sd = nc_model_strong$sigma[2],
                lower.tail = FALSE
            ) # null hypothesis: nc_model_strong > Score$strong
        # lethality
        Score$pvalue_sensitive_lethality <- nc_model_lethality$lambda[1] * pnorm(
            Score$sensitive_lethality,
            mean = nc_model_lethality$mu[1],
            sd = nc_model_lethality$sigma[1],
            lower.tail = FALSE
        ) +
            nc_model_lethality$lambda[2] * pnorm(
                Score$sensitive_lethality,
                mean = nc_model_lethality$mu[2],
                sd = nc_model_lethality$sigma[2],
                lower.tail = FALSE
            ) # null hypothesis: nc_model_lethality > Score$sensitive_lethality
        # recovery
        Score$pvalue_sensitive_recovery <- nc_model_recovery$lambda[1] * pnorm(
            Score$sensitive_recovery,
            mean = nc_model_recovery$mu[1],
            sd = nc_model_recovery$sigma[1],
            lower.tail = FALSE
        ) +
            nc_model_recovery$lambda[2] * pnorm(
                Score$sensitive_recovery,
                mean = nc_model_recovery$mu[2],
                sd = nc_model_recovery$sigma[2],
                lower.tail = FALSE
            ) # null hypothesis: nc_model_recovery > Score$sensitive_recovery
        
        Score$nc_model_strong <- nc_model_strong
        Score$nc_model_lethality <- nc_model_lethality
        Score$nc_model_recovery <- nc_model_recovery
    }
    
    # calculate the adjusted p-values
    if (sum(!is.na(nc_pairs)) != 0) {
        # strong
        Score$fdr_strong <- apply(Score$pvalue_strong, 2, function(x)
            p.adjust(x, method = "fdr"))
        # lethality
        Score$fdr_sensitive_lethality <- apply(Score$pvalue_sensitive_lethality, 2, function(x)
            p.adjust(x, method = "fdr"))
        # recovery
        Score$fdr_sensitive_recovery <- apply(Score$pvalue_sensitive_recovery, 2, function(x)
            p.adjust(x, method = "fdr"))
    }
    
    # output
    Score <- Score[!unlist(lapply(Score, function(x)
        all(is.na(x))))]
    
    return(Score)
}
