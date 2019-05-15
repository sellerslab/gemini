#' Calculate log-fold change
#' @description Given a gemini.input object, calculates log-fold change from \code{counts}.
#' @param Input a gemini.input object containing an object named \code{counts}.
#' @param counts a character indicating the name of a matrix of counts within \code{Input} that can be used to calculate log-fold changes (defaults to "counts").
#' @param sample.column.name a character or integer indicating which column of \code{Input$replicate.map} describes the samples.
#' @param normalize a logical indicating if counts should be median-normalized, see Details (default = TRUE)
#' @param CONSTANT a numeric indicating a constant value that shifts counts to reduce outliers, see Details (default = 32). 
#' @return a gemini object identical to \code{Input} that also contains new objects called \code{LFC} and \code{sample.annot}.
#' 
#' @details 
#' See Methods from Zamanighomi et al. 2019 for a comprehensive 
#' description of the calculation of log-fold change, normalization, 
#' and count processing. 
#'
#' If multiple early time-points are provided for a given sample, they are treated as replicates and averaged,
#' and used to compute log-fold change against any specified late time points.
#' 
#' If a sample has a specific early-time point, these are matched as long as the sample names are identical between
#' the early and late timepoints in \code{sample.column.name}
#'
#' @import magrittr
#' @import dplyr
#' @export
#' 
#' @usage gemini_calculate_lfc(Input, counts = "counts", sample.column.name = "samplename",
#' normalize = TRUE, CONSTANT = 32)
#' 
#' @examples
#' 
#' data("Input", package = "gemini")
#' Input <- gemini_calculate_lfc(Input)
#' 
#' \dontrun{
#' head(Input$LFC)
#' }
#'
gemini_calculate_lfc <- function(Input, 
								 counts = "counts", 
								 sample.column.name = "samplename",
								 normalize = TRUE,
								 CONSTANT = 32){
	stopifnot("gemini.input" %in% class(Input))
	stopifnot(counts %in% names(Input))

	# normalize and scale counts
	mat <- Input[[counts]]
	total.counts = sum(mat, na.rm = TRUE)
	SCALE = (total.counts / ncol(mat))
	if(normalize){
		norm.mat = apply(mat, 2, function(x){
			x = ((x/sum(x, na.rm = TRUE)) * SCALE) + CONSTANT
		})
		Input[[paste0("normalized_",counts)]] <- norm.mat
		
		# compute median-normalized log counts
		data <- norm.mat
	}else{
		data <- mat
	}
	
	data = apply(data, 2, function(x){
		log2(x) - median(log2(x), na.rm = TRUE)
	})

	# compute log-fold changes
	ETP.cols <- which(Input$replicate.map$TP=="ETP")
	ETP.samples <- unique(Input$replicate.map[ETP.cols,][[sample.column.name]])
	if(length(ETP.cols) == 0){
		stop("No ETP samples identified.  Make sure at least one ETP column is specified in Input$replicate.map$TP")
	}else if(length(ETP.cols) == 1){ # If only 1 ETP column specified
		ETP = data[,ETP.cols]
	}else if(length(ETP.cols) > 1 & !any(ETP.samples %in% Input$replicate.map[-ETP.cols,][[sample.column.name]])){ 
	# If multiple ETP replicates belonging to only 1 sample (which is not found in LTP samples, i.e. pDNA)
		ETP = data[,ETP.cols] %>%
			as.data.frame() %>%
			rowMeans(na.rm = TRUE)
	}
	
	# If ETPs match LTPs by sample, that is handled here:
	LTP.cols <- which(Input$replicate.map$TP=="LTP")
	LTP <- as.matrix(data[, LTP.cols])
	colnames(LTP) <- Input$replicate.map$colname[LTP.cols]
	LTP_df <- lapply(unique(Input$replicate.map[[sample.column.name]][Input$replicate.map$TP == "LTP"]), function(x){
		if(!exists('ETP')){ # Check if an ETP dataframe has been established
			ETP.cols = which(Input$replicate.map$TP=="ETP" & Input$replicate.map[[sample.column.name]] == x)
			if(!length(ETP.cols)>0)
				stop("No ETP specified for ", x)
			ETP = data[,ETP.cols] %>%
				as.data.frame() %>%
				rowMeans(na.rm = TRUE) # Create ETP df for this sample
		}
		cols <- Input$replicate.map$colname[Input$replicate.map[[sample.column.name]]==x & Input$replicate.map$TP == "LTP"]
		LFC <- LTP[,match(cols, colnames(LTP), nomatch = 0)] %>%
			as.data.frame(optional = TRUE) %>%
			rowMeans(na.rm = TRUE) %>%
			magrittr::subtract(ETP)
		return(LFC)
	}) %>%
		magrittr::set_names(unique(Input$replicate.map[[sample.column.name]][Input$replicate.map$TP == "LTP"])) %>%
		dplyr::bind_cols() %>%
		as.data.frame(optional = TRUE, stringsAsFactors = FALSE) %>%
		magrittr::set_rownames(Input$guide.pair.annot[,1])

	unique.to.sample <- names(which(apply(Input$replicate.map, 2, function(x) length(unique(x))) == length(unique(Input$replicate.map[[sample.column.name]]))))
	Input$sample.annot <- Input$replicate.map %>%
		dplyr::filter(.$`TP` == "LTP") %>%
		dplyr::select(sample.column.name, 'TP', unique.to.sample) %>%
		unique() %>%
		dplyr::mutate(rowname = colnames(LTP_df))


	Input[["LFC"]] <- LTP_df %>%
		dplyr::select(unique(Input$replicate.map[[sample.column.name]][Input$replicate.map$TP != "ETP"])) %>%
		as.matrix()

	# Return object
	Output <- Input
	class(Output) <- c(class(Output), "gemini.input")
	return(Output)
}
