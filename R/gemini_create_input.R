#' gemini_create_input
#'
#' @description Creates a gemini.input object from a counts matrix with given annotations.
#' @param counts.matrix a matrix of counts with rownames corresponding to features (e.g. guides) and colnames corresponding to samples.
#' @param sample.replicate.annotation a data.frame of annotations for each sample/replicate pair.
#' Note that at least one column in \code{sample.replicate.annotation} must correspond to the colnames of \code{counts.matrix} (see Details) (default = NULL)
#' @param guide.annotation a data.frame of annotations for each guide.  Note that at least one column in \code{guide.annotation} must correspond to the rownames of counts.matrix (default = NULL)
#' @param samplesAreColumns a logical indicating if samples are on the columns or rows of counts.matrix. (default = TRUE)
#' @param sample.column.name a character or integer indicating which column of \code{sample.replicate.annotation} describes the samples.
#' @param gene.column.names a character or integer vector of length(2) indicating which columns of \code{guide.annotation} describe the genes being targeted.
#' @param ETP.column a character or integer vector indicating which column(s) of \code{counts.matrix} contain the early time-point(s) of the screen (i.e. pDNA, early sequencing, etc.).  Defaults to the first column.
#' @param LTP.column a character or integer vector indicating which column(s) is the later time-point of the screen (i.e. day21, post-treatment, etc.).  Defaults to \code{(1:ncol(counts.matrix))[-ETP.column]}, or all other columns except for those specified by \code{ETP.column}.
#' @param verbose Verbosity (default FALSE)
#' @return a gemini.input object
#'
#' @details
#' This function initializes a gemini.input object from a counts matrix. There are a few key assumptions made in the input format.
#' \itemize{
#' \item The counts matrix is regular.
#' \item The counts matrix structure is in accordance with the \code{samplesAreColumns} parameter.
#' \item The first column of \code{sample.replicate.annotation} matches with the existing dimension names of the counts matrix.
#' \item The first column of \code{guide.annotations} matches with the existing dimension names of the counts matrix.
#' \item \code{sample.column.name} must specify a column in \code{sample.replicate.annotation} (either by name or index) that describes unique samples.
#' \item \code{gene.column.names} must specify two columns in \code{sample.replicate.annotation} (either by name or index) that describe genes.
#' }
#'
#' @importFrom dplyr mutate
#'
#' @export
gemini_create_input <-
	function(counts.matrix,
			 sample.replicate.annotation = NULL,
			 guide.annotation = NULL,
			 samplesAreColumns = TRUE,
			 sample.column.name = NULL,
			 gene.column.names = NULL,
			 ETP.column = 1,
			 LTP.column = NULL,
			 verbose = F) {
		
	  # Check ETP/LTP column identification
	  if(is.numeric(ETP.column) & is.null(LTP.column)){
	    LTP.column <- (1:ncol(counts.matrix))[-ETP.column]
	  }else if(is.character(ETP.column) & is.null(LTP.column)){
	    ETP.column <- which(colnames(counts.matrix) %in% ETP.column)
	    LTP.column <- (1:ncol(counts.matrix))[-ETP.column]
	  }
	  
		# Require dimension names for counts matrix if no guide and replicate annotations provided
		if (is.null(dimnames(counts.matrix)) | is.null(guide.annotation) | is.null(sample.replicate.annotation))
			stop("No dimnames for counts.matrix - no annotations available.", "")
	  
	  # Require sample.column.name and gene.column.names specification
	  if (is.null(gene.column.names) | is.null(sample.column.name)){
	    stop("Did you provide gene.column.names and/or sample.column.name?")
	  }
		
		# transpose matrix
		if (!samplesAreColumns) {
		  if(verbose) message("Transposing matrix...")
			# transpose and preserve dimnames
			dn <- dimnames(counts.matrix)
			counts.matrix <- t(counts.matrix)
			dimnames(counts.matrix) <- rev(dn)
		}
		
		# default guide annotations to rownames of counts matrix
		gannot <-
			data.frame(rowname = rownames(counts.matrix),
					   stringsAsFactors = F)
		
		# Default sample annotations to column names of counts matrix, ordering by ETP -> LTP
		sannot <-
			data.frame(colname = colnames(counts.matrix)[c(ETP.column, LTP.column)],
					   stringsAsFactors = F, row.names = 1:length(c(ETP.column, LTP.column))) %>%
			dplyr::mutate(TP = c(rep("ETP", length(ETP.column)), rep("LTP", length(LTP.column))))
		
		# Merge existing sample annotations with colnames, ensuring formatting and matching names
		if (!is.null(sample.replicate.annotation) & !is.null(sample.column.name)) {
			colnames(sample.replicate.annotation)[colnames(sample.replicate.annotation) == sample.column.name] <- "samplename" # Set sample column name to "samplename"
			if(verbose) message("Merging sample annotations with colnames of counts.matrix...")
			i = which(apply(sample.replicate.annotation, 2, function(x)
				all(x %in% sannot[, 1])))
			if(!length(i) > 0){
			  if(verbose) message("No columns found in sample.replicate.annotation which completely match colnames of counts.matrix...")
			}
			sannot <- merge(sannot, sample.replicate.annotation, by.x = 1, by.y = i[1], no.dups = F, all = F, sort = F, suffixes = c("", ".y"))
		}else{
			stop("Could not determine samplename.  Please add sample/replicate annotation and specify and sample.column.name.  See ?gemini_create_input.")
		}
		
		# Merge guide annotations with existing rownames, ensuring formatting and matching names
		if (!is.null(guide.annotation)) {
		  if(verbose) message("Merging guide annotations with rownames()...")
			i = which(apply(guide.annotation, 2, function(x)
				all(x %in% gannot[, 1])))
			if(!length(i) > 0){
			  if(verbose) message("No columns found in guide.annotation which completely match rownames()...")
			}
			gannot <- merge(gannot, guide.annotation, by.x = 1, by.y = i, no.dups = F, all = F, sort = F, suffixes = c("", ".y"))
		}else{
			stop("Could not determine gene/guide data.  Please add guide annotation and specify and gene.column.names. See ?gemini_create_input.")
		}
		
		# Create new Input object
		Output <- list(
			counts = data.matrix(counts.matrix[, c(ETP.column, LTP.column)]),
			replicate.map = as.data.frame(sannot, optional = T, row.names = 1:nrow(sannot)),
			guide.pair.annot = as.data.frame(gannot, optional = T, rownames = 1:nrow(gannot))
		)
		
		Output <- gemini_prepare_input(Output, gene.columns = gene.column.names)
		
		class(Output) <- c(class(Output), "gemini.input")
		if(verbose) message("Created gemini input object.")
		return(Output)
	}
