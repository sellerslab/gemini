#' Big Papi counts matrix
#' @name counts
#' @description Counts matrix from the Big Papi SynLet dataset. 
#' The Big Papi dataset was published in Najm et al. 2018 (doi: 10.1038/nbt.4048). 
#' The counts were acquired from Supplementary Table 3.
#' @docType data
#' @seealso https://www.nature.com/articles/nbt.4048
#' @keywords data
"counts"

#' Big Papi guide annotations
#'
#' @name guide.annotation
#' 
#' @description Guide and gene annotations for the Big Papi SynLet dataset. 
#' The Big Papi dataset was published in Najm et al. 2018 (doi: 10.1038/nbt.4048). 
#' The guide and gene annotations were acquired from Supplementary Table 2.
#' @docType data
#' @seealso https://www.nature.com/articles/nbt.4048
#' @keywords data
"guide.annotation"

#' Big Papi sample and replicate annotations
#'
#' @name sample.replicate.annotation
#' @description Sample and replicate annotations for the Big Papi SynLet dataset.
#' The Big Papi dataset was published in Najm et al. 2018 (doi: 10.1038/nbt.4048). 
#' The sample and replicate annotations were acquired from Supplementary Table 3
#' @docType data
#' @seealso https://www.nature.com/articles/nbt.4048
#' @keywords data
"sample.replicate.annotation"

#' Input object from Big Papi
#'
#' @name Input
#' @description A gemini.input object created from the Big Papi SynLet dataset.
#' The Big Papi dataset was published in Najm et al. 2018 (doi: 10.1038/nbt.4048).
#' The Input object here is created from the data in \link{counts}, 
#' \link{guide.annotation}, and \link{sample.replicate.annotation} 
#' using the \code{\link{gemini_create_input}} function.
#' @docType data
#' @seealso https://www.nature.com/articles/nbt.4048
#' @keywords data
"Input"

#' Model object from Big Papi
#'
#' @name Model
#' @description A gemini.model object created from the Big Papi SynLet dataset after \code{gemini_inference} was run.
#' The Big Papi dataset was published in Najm et al. 2018 (doi: 10.1038/nbt.4048).
#' The Model object here is created from the data in \link{Input} using the \code{\link{gemini_initialize}} function.
#' @docType data
#' @seealso https://www.nature.com/articles/nbt.4048
#' @keywords data
"Model"
