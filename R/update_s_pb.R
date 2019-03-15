#' update_s
#' @description Update values of s using data from \code{Input} and current values of other parameters.
#'
#' @param Model a Model object of class gemini.model
#' @param Input an Input object of class gemini
#' @param LFC an object within \code{Input} containing log-fold change values
#' @param guide.pair.annot an object within \code{Input} mapping guide pairs to individual genes.
#' @param mean_s numeric indicating prior mean of s (default 0)
#' @param sd_s numeric indicating prior sd of s (default 10)
#' @param cores a numeric indicating the number of cores to use, see \code{\link{gemini_parallel}}. default=1.
#' @param verbose default FALSE
#'
#' @return An object of class gemini.model
#'
#' @import pbmcapply
#' @export
#'
update_s_pb <- function(Model, 
                        Input, 
                        LFC = "LFC", 
                        guide.pair.annot = "guide.pair.annot", 
                        mean_s = 0, 
                        sd_s = 10, 
                        cores = 1, 
                        verbose = F){
  
  if(verbose) {
    message("Updating s...")
    message("\tUsing ", cores, " core(s).")
    tstart = Sys.time()
  }
  
  LFC <- Input[[LFC]]
  guide2gene <- Input[[guide.pair.annot]]
  
  # mean of gamma distribution
  tau = Model$alpha/Model$beta
  
  # updating synergy
  s_loop <- function(i, Model, guide2gene, LFC, tau, mean_s, sd_s){
    gihj = Model$hash_s[[i]]
    gi = Model$hashes_x$paired_guide[gihj,1]
    hj = Model$hashes_x$paired_guide[gihj,2]
    g = guide2gene[match(gihj,guide2gene[,1]),2]
    h = guide2gene[match(gihj,guide2gene[,1]),3]
    
    # numerator for pair of genes
    numerator = (Model$xx[gihj]*tau[gihj,])*(LFC[gihj,] - Model$x_seq1[gi]*Model$y[g,] - Model$x_seq2[hj]*Model$y[h,])
    numerator = mean_s/(sd_s^2) + colSums(as.matrix(numerator), na.rm = T)
    
    # denominator for pair of genes
    denominator = 1/(sd_s^2) + colSums(as.matrix(Model$xx2[gihj]*tau[gihj,]), na.rm = T)
    
    # updating s and s2
    s = numerator/denominator
    s2 = s^2 + 1/denominator
    
    return(list(s = s, s2 = s2))
  }
  
  s_res <- pbmcapply::pbmclapply(X = rownames(Model$s), FUN = s_loop, Model = Model, 
                    guide2gene = guide2gene, LFC = LFC, tau = tau,
                    mean_s = mean_s, sd_s = sd_s, mc.cores = cores, max.vector.size = 9999)
  
  Model$s[] <- lapply(s_res, extract, "s") %>%
    unlist(recursive = F, use.names = F) %>%
    do.call(rbind, .)
  
  Model$s2[] <- lapply(s_res, extract, "s2") %>%
    unlist(recursive = F, use.names = F) %>%
    do.call(rbind, .)

  # output
  if(verbose){
    tend = Sys.time()
    tdiff = difftime(tend, tstart)
    message("\tCompleted update of s.")
    message("\tTime to completion: ", round(tdiff, digits = 3), ' ', units(tdiff))
  }
  return(Model)
}
