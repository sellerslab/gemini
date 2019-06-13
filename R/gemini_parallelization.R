#' @title Parallelization in GEMINI
#' @name gemini_parallelization
#' @description Notes about parallelization in combinatorial CRISPR analysis
#' 
#' @section Implementation:
#' 
#' To improve efficiency and scalability, GEMINI employs parallelization, enabling a rapid initialization and update routine. 
#' Parallelization was implemented using the R \code{pbmclapply} package, and specifically through the \code{\link[pbmcapply]{pbmclapply}} function. 
#' As pbmclapply (and it's parent function, \code{\link[parallel]{mclapply}}) relies on forking, parallelization is currently limited to Unix-like (Linux flavors and macOS) machines, but may be extended to other OS in later versions through socket parallelization. In the event that GEMINI is used on a Windows machine, all computations are performed in serial.
#' 
#' @section Parallelized processes:
#' 
#' GEMINI enables parallelization of both initialization and inference. For initialization, parallelization is used to quickly hash the data. For the inference procedure, parallelization is used to speed up the CAVI approach by performing updates independently across cores. 
#' 
#' @section Caveats:
#' To note, there is usually a trade-off in terms of number of cores and shared resources/memory transfer. See inst/figs/gemini-benchmarking.png for a visual depiction.
#' 
#' Also, while most functions in GEMINI have been parallelized, initialization of tau, update of x (group 3), and updates of y are performed in serial. As such, especially in large libraries, tau initialization and x (group 3) updates may take some time. Updating y is usually fast, as the number of genes tends to be the smallest in parameter space. However, these functions may be parallelized in the future as well.
#' 
#' @section Active Development:
#' Noticed that in some cases, after using multicore processing through pbmclapply,
#' warnings have been produced: 
#' "In selectChildren(pids[!fin], -1) : cannot wait for child ... as it does not exist"
#' "In parallel::mccollect(...) : 1 parallel job did not deliver a result"
#' R.version=3.6.0
#' These warnings, although they appear menacing, are in fact harmless. 
#' See: http://r.789695.n4.nabble.com/Strange-error-messages-from-parallel-mcparallel-family-under-3-6-0-td4756875.html#a4756939
#' 
#' 
NULL
