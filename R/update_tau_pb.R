#' update_tau_pb
#' @description Update parameters of tau using data from \code{Input} and current values of other parameters.
#'
#' @param Model a Model object of class gemini.model
#' @param Input an Input object of class gemini.input
#' @param LFC an object within \code{Input} containing log-fold change values
#' @param guide.pair.annot an object within \code{Input} mapping guide pairs to individual genes.
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
update_tau_pb <- function(Model,
                          Input,
                          LFC = "LFC",
                          guide.pair.annot = "guide.pair.annot",
                          nc_gene,
                          pattern_join = ";",
                          cores = 1,
                          verbose = F){
  if(verbose){
    message("Updating tau...")
    message("\tUsing ", cores, ' core(s).')
    tstart = Sys.time()
  } 
	
  LFC <- Input[[LFC]]
  guide2gene <- Input[[guide.pair.annot]]
  
  # mean of gamma distribution
  tau = Model$alpha/Model$beta
  
  # define loop function
  tau_loop <- function(i, Model, guide2gene, LFC, pattern_join, nc_gene){
    gi = Model$hashes_x$paired_guide[i,1]
    hj = Model$hashes_x$paired_guide[i,2]
    g = guide2gene[match(i,guide2gene[,1]),2]
    h = guide2gene[match(i,guide2gene[,1]),3]
    gh = paste(sort(c(g,h)),collapse = pattern_join)

    # calculating beta*
    gnc = g %in% nc_gene
    hnc = h %in% nc_gene
    if (!(gnc) & !(hnc)){
      beta_star = LFC[i,]^2 - 2*LFC[i,]*(Model$x_seq1[gi]*Model$y[g,] + Model$x_seq2[hj]*Model$y[h,] + Model$xx[i]*Model$s[gh,]) + Model$x2_seq1[gi]*Model$y2[g,] + 2*Model$x_seq1[gi]*Model$y[g,]*(Model$x_seq2[hj]*Model$y[h,] + Model$xx[i]*Model$s[gh,]) + Model$x2_seq2[hj]*Model$y2[h,] + Model$xx2[i]*Model$s2[gh,] + 2*Model$x_seq2[hj]*Model$y[h,]*Model$xx[i]*Model$s[gh,]
    } else if (gnc & !(hnc)){
      beta_star = LFC[i,]^2 + Model$x2_seq2[hj]*Model$y2[h,] - 2*LFC[i,]*Model$x_seq2[hj]*Model$y[h,]
    } else if (!(gnc) & hnc){
      beta_star = LFC[i,]^2 + Model$x2_seq1[gi]*Model$y2[g,] - 2*LFC[i,]*Model$x_seq1[gi]*Model$y[g,]
    } else {
      beta_star = LFC[i,]^2
    }

    # updating alpha and beta
    alpha_i = Model$prior_shape + 0.5
    beta_i = Model$beta_prior[i,] + 0.5*beta_star
    return(list(alpha = alpha_i, beta = beta_i))
    
  }
  
  res <- pbmcapply::pbmclapply(X = rownames(Model$alpha), FUN = tau_loop, 
                           Model = Model, guide2gene = guide2gene, LFC = LFC, 
                           pattern_join = pattern_join, nc_gene = nc_gene, mc.cores = cores, max.vector.size = 9999)
  
  Model$beta[,] <- lapply(res, magrittr::extract, "beta") %>%
    unlist(recursive = F, use.names = F) %>%
    do.call(rbind, .)
  
  Model$alpha[,] <- lapply(res, magrittr::extract, "alpha") %>%
  	unlist(recursive = F, use.names = F) %>%
  	do.call(rbind, .)
  
  # output
  if(verbose){
    tend = Sys.time()
    tdiff = difftime(tend, tstart)
    message("\tCompleted update of tau.")
    message("\tTime to completion: ", round(tdiff, digits = 3), ' ', units(tdiff))
  }
  return(Model)
}


