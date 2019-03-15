#' @title Parallelization in GEMINI
#' @name gemini_parallel
#' 
#' @section Implementation:
#' 
#' To improve efficiency and scalability, GEMINI employs parallelization, enabling a rapid initialization and update routine. 
#' This parallelization was implemented using the R \code{pbmclapply} package, and specifically through the \code{\link[pbmcapply]{pbmclapply}} function. 
#' As pbmclapply (and it's parent function, \code{\link[parallel]{mclapply}}) relies on forking, parallelization is currently limited to Unix-like (Linux flavors and macOS) machines, but may be extended to other OS in later versions through socket parallelization.
#' 
#' @section Caveats:
#' To note, there is usually a trade-off in terms of number of cores and shared resources/memory transfer. Please see the README.html for a graphical depiction of this trade-off.
#' 
#' Finally, while most functions in GEMINI have been parallelized, initialization of tau and updates of y are not parallelizable.  As such, especially in large libraries, tau initialization may take some time. Updating y is usually fast, as the number of genes tends to be the smallest in parameter space. However, these functions may be parallelized in the future as well.
#' 
NULL