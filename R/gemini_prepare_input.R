#' Prepare input before Model creation
#'
#' @description This is an internal function to GEMINI, allowing for data cleanup and preprocessing before a Model object is created. This removes any gene pairs targeting the same gene twice, and removes any empty samples/replicates.
#'
#' @param Input An object of class gemini.input
#' @param gene.columns a character vector of length(2)
#' @param sample.col.name a character indicating the name of the sample column (default = "samplename")
#'
#' @return a (prepared) gemini.input object
#' 
#' @examples 
#' Input %<>% gemini_prepare_input(gene.columns = c("U6.gene", "H1.gene"))
#' 
#' @export
#'

gemini_prepare_input <-
	function(Input, gene.columns, sample.col.name = "samplename") {
		LFC = "LFC"
		newInput <- list()
		class(newInput) <- base::union(class(Input), "gemini.input")
		
		dups <-
			apply(Input$guide.pair.annot[, gene.columns], 1, function(x)
				any(duplicated(x)))
		newInput[["counts"]] <-
			Input$counts[!dups, colSums(!is.na(Input$counts)) != 0]
		newInput[["replicate.map"]] <- Input$replicate.map %>%
			dplyr::filter(colSums(!is.na(Input$counts)) != 0) %>%
			dplyr::select(c("colname", sample.col.name, "TP"))
		newInput[["guide.pair.annot"]] <- Input$guide.pair.annot %>%
			dplyr::filter(!dups) %>%
			dplyr::select(c("rowname", gene.columns))
		if (LFC %in% names(Input)) {
			newInput[["LFC"]] <-
				data.matrix(Input[[LFC]][!dups, unique(newInput$replicate.map$sample[newInput$replicate.map$TP ==
																					 	"LTP"])])
			colnames(newInput[["LFC"]]) <- colnames(Input[[LFC]])
		}
		return(newInput)
	}
