#' gemini_boxplot
#'
#' @description A function to visualize the results of GEMINI over the raw data.
#'
#' @param Model a gemini.model object
#' @param g a character naming a gene to visualize
#' @param h a character naming another gene to visualize
#' @param sample a character naming the sample to visualize
#' @param show_inference a logical indicating whether to show the 
#' inferred individual or combined values for each gene/gene pair (default TRUE)
#' @param color_x a logical indicating whether to visualize the 
#' sample-independent effects for each individual guide or guide pair (default FALSE)
#' @param identify_guides a logical indicating whether to identify 
#' guides with unique colors and shapes (default FALSE)
#' @param nc_gene a character naming the gene to use as a negative control, 
#' to be paired with each individual g and h. Defaults to \code{Model$nc_gene}.
#'
#' @details Raw LFC data is plotted for each gene combination (`g`-`nc_gene`, `h`-`nc_gene`, `g`-`h`) in a standard boxplot.
#' Horizontal green line segments are plotted over the box plots indicating the individual gene effects or
#' the inferred total effect of a particular gene combination. Each guide
#' pair can be colored based on the inferred sample independent effects
#' \eqn{g_i}, \eqn{h_j}, and \eqn{g_i,h_j}. Additionally, colors and shapes 
#' can be used to distinguish unique guides targeting gene g and h, respectively.
#'
#' @import ggplot2
#' @importFrom dplyr mutate
#' @importFrom grDevices hcl
#'
#' @return a ggplot2 object
#'
#' @export
#'
#' @examples
#'
#' data("Model", package = "gemini")
#'
#' gemini_boxplot(Model, g = "BRCA2", h = "PARP1", nc_gene = "CD81", 
#' sample = "A549", show_inference = TRUE,
#' color_x = FALSE, identify_guides = FALSE)
#'
#'
gemini_boxplot <- function(Model,
                           g,
                           h,
                           nc_gene = NULL,
                           sample,
                           show_inference = TRUE,
                           color_x = FALSE,
                           identify_guides = FALSE) {
    # Extract params from Model
    Input = Model$Input
    if (is.null(nc_gene))
        nc_gene = Model$nc_gene
    stopifnot(!is.null(nc_gene)) # if no nc_gene provided, require for plotting.
    
    # load required packages, and suggest installation if not already installed.
    if (!requireNamespace("ggplot2"))
        stop("Please install the ggplot2 package: install.packages('ggplot2')")
    if (!requireNamespace("magrittr"))
        stop("Please install the magrittr package: install.packages('magrittr')")
    
    # Check inputs:
    stopifnot(sample %in% colnames(Model$y))
    stopifnot(g %in% rownames(Model$y) & h %in% rownames(Model$y))
    if (color_x & identify_guides) {
        color_x = FALSE
        warning(
            "color_x and identify_guides cannot be concurrently visualized. Defaulting to identify_guides."
        )
    }
    
    gihj = Model$hash_s[[paste0(sort(c(g, h)), collapse = Model$pattern_join)]]
    gh = paste0(sort(c(g, h)), collapse = Model$pattern_join)
    ginc = Model$hash_s[[paste0(sort(c(g, nc_gene)), collapse = Model$pattern_join)]]
    hjnc = Model$hash_s[[paste0(sort(c(h, nc_gene)), collapse = Model$pattern_join)]]
    
    # For cases in which the negative control pairs were not 
    # included in the inference (i.e. used in the inference process)
    if (is.null(ginc) | is.null(hjnc)) {
        ginc = Input$guide.pair.annot$rowname[Input$guide.pair.annot[, 2] == g &
                                                  Input$guide.pair.annot[, 3] == nc_gene |
                                                  Input$guide.pair.annot[, 2] == nc_gene &
                                                  Input$guide.pair.annot[, 3] == g]
        hjnc = Input$guide.pair.annot$rowname[Input$guide.pair.annot[, 2] == h &
                                                  Input$guide.pair.annot[, 3] == nc_gene |
                                                  Input$guide.pair.annot[, 2] == nc_gene &
                                                  Input$guide.pair.annot[, 3] == h]
    }
    
    # Create data frame using negative control pairs with all g_i
    df_ginc <- lapply(ginc, function(x) {
        df = data.frame(
            gihj = x,
            gi = strsplit(x, split = Model$pattern_split, fixed = TRUE)[[1]][1],
            hj = strsplit(x, split = Model$pattern_split, fixed = TRUE)[[1]][2],
            D = Input$LFC[x, sample],
            stringsAsFactors = FALSE
        ) %>%
            mutate(x_gi = Model$x[.$`gi`]) %>%
            mutate(x_hj = Model$x[.$`hj`]) %>%
            mutate(xx_gihj = Model$xx[x]) %>%
            mutate(y = Model$y[g, sample]) %>%
            mutate(label = paste0(c(g, nc_gene), collapse = Model$pattern_join))
    }) %>%
        do.call(rbind, .)
    
    # Create data frame using negative control pairs with all h_j
    df_hjnc <- lapply(hjnc, function(x) {
        df = data.frame(
            gihj = x,
            gi = strsplit(x, split = Model$pattern_split, fixed = TRUE)[[1]][1],
            hj = strsplit(x, split = Model$pattern_split, fixed = TRUE)[[1]][2],
            stringsAsFactors = FALSE
        ) %>%
            mutate(x_gi = Model$x[gi]) %>%
            mutate(x_hj = Model$x[hj]) %>%
            mutate(xx_gihj = Model$xx[x]) %>%
            mutate(D = Input$LFC[x, sample]) %>%
            mutate(y = Model$y[h, sample]) %>%
            mutate(label = paste0(c(h, nc_gene), collapse = Model$pattern_join))
    }) %>%
        do.call(rbind, .)
    
    # Create data frame using all g_i with all h_j
    df_gihj <- lapply(gihj, function(x) {
        df = data.frame(
            gihj = x,
            gi = strsplit(x, split = Model$pattern_split, fixed = TRUE)[[1]][1],
            hj = strsplit(x, split = Model$pattern_split, fixed = TRUE)[[1]][2],
            stringsAsFactors = FALSE
        ) %>%
            mutate(x_gi = Model$x[.$`gi`]) %>%
            mutate(x_hj = Model$x[.$`hj`]) %>%
            mutate(xx_gihj = Model$xx[x]) %>%
            mutate(D = Input$LFC[x, sample]) %>%
            mutate(y = Model$y[g, sample] + Model$y[h, sample] + Model$s[gh, sample]) %>%
            mutate(label = gh)
    }) %>%
        do.call(rbind, .)
    
    # Bind all dataframes together to make final data structure for plotting
    data <- do.call(rbind, list(df_ginc, df_hjnc, df_gihj))
    
    # Add labels to x-axis of boxplot
    data$label <-
        factor(data$label, levels = c(
            paste0(c(g, nc_gene), collapse = Model$pattern_join),
            gh,
            paste0(c(h, nc_gene), collapse = Model$pattern_join)
        ))
    
    if (color_x) {
        # Process x-values for visualization
        xs = cbind(data$x_gi, data$x_hj, data$xx_gihj)
        abs_xs = (xs - 1)
        
        # data %<>%
        #     mutate(not_NA = apply(xs, 1, function(x)
        #         sum(!is.na(x)))) %>%
        #     mutate(abs_max = apply(abs_xs, 1, function(x)
        #         x[which.max(abs(x))])) %>%
        #     mutate(avg_x = apply(xs, 1, mean, na.rm = TRUE))
        
        data$vis_x <- data$xx_gihj
        data$vis_x[is.na(data$vis_x)] <- data$x_gi[is.na(data$vis_x)]
        data$vis_x[is.na(data$vis_x)] <- data$x_hj[is.na(data$vis_x)]
        
        # Generate plot using data
        p = ggplot(data = data, aes_string(
            x = "label",
            y = "D",
            color = "vis_x"
        )) +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter(
                data = data,
                mapping = aes_string(color = "vis_x"),
                width = 0.2,
                size = 2.5,
                height = 0
            ) +
            scale_color_gradient2(
                midpoint = 1,
                high = "red",
                mid = "black",
                low = 'yellow'
            ) +
            labs(
                x = "",
                y = "Log-fold change",
                color = "Avg inferred x",
                title = sample
            ) +
            theme(
                axis.text = element_text(color = 'black'),
                axis.ticks = element_line(color = 'black'),
                panel.background = element_blank()
            ) +
            labs(x = "",
                 y = "Log-fold change",
                 title = sample)
    } else if (identify_guides) {
        # Get guide sequences
        g_seqs <-
            names(which(Model$hashes_x$gene_hash == g))
        h_seqs <-
            names(which(Model$hashes_x$gene_hash == h))
        # Identify color/shape matches for each guide
        data$color <- g_seqs[match(data$gi, g_seqs)]
        data$color[is.na(data$color)] <-
            g_seqs[match(data$hj[is.na(data$color)], g_seqs)]
        data$color[is.na(data$color)] <- nc_gene
        
        data$shape <- h_seqs[match(data$gi, h_seqs)]
        data$shape[is.na(data$shape)] <-
            h_seqs[match(data$hj[is.na(data$shape)], h_seqs)]
        data$shape[is.na(data$shape)] <- nc_gene
        
        # Transform to factor for plotting
        data$shape <-
            factor(data$shape, levels = c(unique(data$shape[data$shape != nc_gene]), nc_gene))
        nshapes = length(unique(data$shape[!is.na(data$shape)])) - 1
        
        data$color <-
            factor(data$color, levels = c(unique(data$color[data$color != nc_gene]), nc_gene))
        ncolors = length(unique(data$color[!is.na(data$color)])) - 1
        
        # Create color picking function to select from color wheel
        gg_color_hue <- function(n) {
            hues = seq(15, 375, length = n + 1)
            pickedcolors = grDevices::hcl(h = hues,
                                          l = 65,
                                          c = 100)[seq_len(n)]
            pickedcolors = c(pickedcolors, grDevices::hcl(0, 0, 0))
            return(pickedcolors)
        }
        
        gg_shape_select <- function(n) {
            allshapes = c(3, 7, 8, 15, 18, 17, 11, 9, 10, 4)
            if (n > length(allshapes)) {
                warning(
                    "More guides than available shapes: Not all guides (shapes) may be distinguishable!"
                )
                pickedshapes <- sample(allshapes,
                                       size = n,
                                       replace = TRUE)
            } else{
                pickedshapes <- sample(allshapes,
                                       size = n,
                                       replace = FALSE)
            }
            pickedshapes <- c(pickedshapes, 16)
            return(pickedshapes)
        }
        
        p = ggplot(data = data, aes_string(x = "label", y = "D")) +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter(
                data = data,
                aes_string(color = "color", shape = "shape"),
                width = 0.2,
                size = 2.5,
                height = 0
            ) +
            scale_color_manual(values = gg_color_hue(ncolors)) +
            scale_shape_manual(values = gg_shape_select(nshapes)) +
            theme(
                axis.text = element_text(color = 'black'),
                axis.ticks = element_line(color = 'black'),
                panel.background = element_blank()
            ) +
            labs(x = "",
                 y = "Log-fold change",
                 title = sample)
        
        
    } else{
        # Plot default boxplot values
        p = ggplot(data = data, aes_string(x = "label", y = "D")) +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter(
                data = data,
                width = 0.2,
                size = 2.5,
                height = 0,
                color = 'black'
            ) +
            theme(
                axis.text = element_text(color = 'black'),
                axis.ticks = element_line(color = 'black'),
                panel.background = element_blank()
            ) +
            labs(x = "",
                 y = "Log-fold change",
                 title = sample)
    }
    
    if (show_inference) {
        p = p + geom_segment(
            mapping = aes(
                x = as.numeric(label) - 0.5,
                y = y,
                yend = y,
                xend = as.numeric(label) + 0.5
            ),
            color = "green"
        )
    }
    
    return(p)
}
