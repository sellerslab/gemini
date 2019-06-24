#' update_x_p
#' @description Update values of x using data from \code{Input} and current values of other parameters.
#'
#' @param Model a Model object of class gemini.model
#' @param mean_x a numeric indicating prior mean of x
#' @param sd_x a numeric indicating prior sd of x
#' @param mean_xx a numeric indicating prior mean of xx
#' @param sd_xx a numeric indicating prior sd of xx
#' @param cores a numeric indicating the number of cores to use.  See \code{\link[gemini]{gemini_parallelization}} for details.  (default=1).
#' @param verbose default FALSE
#'
#' @return An object of class gemini.model
#' @note The structure of the screen may impede parallelization. Our ability to parallelize the updates is 
#' contingent upon the independence between guides in position 1 and 2.
#' To account for potential dependence, we define three groups of guides: group 1 (guides only in position 1
#' ), group 2 (guides only in position 2), and group 3 (guides in both position 1 and position 2). 
#' Parallelization is possible for groups 1 and 2, but only serial updates are possible for group 3. 
#' As such, updates for this group will take longer.
#' @importFrom parallel mclapply
#' @importFrom pbmcapply pbmclapply
#' @importFrom magrittr extract
#'
#' @importFrom utils setTxtProgressBar
#' @importFrom utils getTxtProgressBar
#'
#' @examples 
#' data("Model", package = "gemini")
#' Model %<>% update_x_pb()
#'
#' @export
update_x_pb <- function(Model,
                        mean_x = 1,
                        sd_x = 1,
                        mean_xx = 1,
                        sd_xx = 1,
                        cores = 1,
                        verbose = FALSE) {
    if (verbose) {
        message("Updating x...")
        message("\tUsing ", cores, " core(s).")
        tstart = Sys.time()
    }
    
    Input <- Model$Input
    guide2gene <- Input[[Model$guide.pair.annot]]
    LFC <- Input[[Model$LFC.name]]
    
    # genes corresponding to x_seq
    x_seq_genes = Model$hashes_x$gene_hash
    
    # mean of gamma distribution
    tau = Model$alpha / Model$beta
    
    # Allow looping through sd_x if lists are provided
    if (length(sd_x) == 1) {
        sd_x = rep(sd_x, length(Model$hashes_x$hash))
        names(sd_x) <- names(Model$hashes_x$hash)
    }
    if (length(mean_x) == 1) {
        mean_x = rep(mean_x, length(Model$hashes_x$hash))
        names(mean_x) <- names(Model$hashes_x$hash)
    }
    if (length(mean_xx) == 1) {
        mean_xx = rep(mean_xx, nrow(Model$hashes_x$paired_guide))
        names(mean_xx) <- rownames(Model$hashes_x$paired_guide)
    }
    if (length(sd_xx) == 1) {
        sd_xx = rep(sd_xx, nrow(Model$hashes_x$paired_guide))
        names(sd_xx) <- rownames(Model$hashes_x$paired_guide)
    }
    
    # update x_seq1 and x2_seq1
    x_loop <-
        function(i,
                 Model,
                 guide2gene,
                 LFC,
                 tau,
                 x_seq_genes,
                 mean_x,
                 sd_x) {
            # get individual mean_xi and sd_xi for guide i
            # NOTE: mean_x/sd_x should be a named numeric vector, where names correspond to guide labels
            mean_xi = mean_x[i]
            sd_xi = sd_x[i]
            # x_seq to gene names
            g = as.character(x_seq_genes[i]) # gene g corresponding to guide i
            gihj = Model$hashes_x$hash[[i]] # all guide pairs containing i
            h.col = vapply(gihj, function(x){
                guidepair = strsplit(x,split = Model$pattern_split, fixed = TRUE)[[1]]
                return(which(guidepair!=i))
            }, numeric(1))
            
            # identify genes h paired with gene g
            h = guide2gene[match(gihj, guide2gene[, 1]), 2] 
            h[h.col==2] <- guide2gene[match(gihj, guide2gene[, 1]), 3][h.col==2] # if any h in position 1
            
            # identify guides hj paired with guide gi
            hj = Model$hashes_x$paired_guide[match(gihj, rownames(Model$hashes_x$paired_guide)), 1] # guides for genes h
            hj[h.col==2] <- Model$hashes_x$paired_guide[match(gihj, rownames(Model$hashes_x$paired_guide)), 2][h.col==2] # if any hj in position 1
            
            gh = apply(cbind(rep(g, length(h)), h), 1, function(x)
                paste(sort(x), collapse = Model$pattern_join)) # gene g paired with genes h
            # calculating numerator
            numerator = LFC[gihj[!h %in% Model$nc_gene], ] - Model$x[hj[!h %in% Model$nc_gene]] *
                Model$y[h[!h %in% Model$nc_gene], ] - Model$xx[gihj[!h %in% Model$nc_gene]] * Model$s[gh[!h %in% Model$nc_gene], ]
            numerator = tau[gihj[!h %in% Model$nc_gene], ] * numerator
            if (sum(!h %in% Model$nc_gene)>1 | sum(!h %in% Model$nc_gene)==0){
                numerator = colSums(as.matrix(numerator), na.rm = TRUE)    
            } else{
                numerator = colSums(as.matrix(t(numerator)), na.rm = TRUE)
            }
            if (sum(h %in% Model$nc_gene)>1 | sum(h %in% Model$nc_gene)==0){
                numerator = numerator +
                    colSums(as.matrix(tau[gihj[h %in% Model$nc_gene], ] * LFC[gihj[h %in% Model$nc_gene], ]), na.rm = TRUE)
            } else{
                numerator = numerator +
                    colSums(as.matrix(t(tau[gihj[h %in% Model$nc_gene], ] * LFC[gihj[h %in% Model$nc_gene], ])), na.rm = TRUE)
            }
            #numerator = colSums(as.matrix(numerator), na.rm = TRUE) +
            #    colSums(as.matrix(tau[gihj[h %in% Model$nc_gene], ] * LFC[gihj[h %in% Model$nc_gene], ]), na.rm = TRUE)
            numerator = sum(Model$y[g, ] * numerator, na.rm = TRUE)
            numerator = mean_xi / (sd_xi ^ 2) + numerator
            # calculating denominator
            if(length(gihj)>1){
                denominator = sum(colSums(as.matrix(tau[gihj, ]), na.rm = TRUE) * Model$y2[g, ], na.rm = TRUE)
            } else{
                denominator = sum(colSums(as.matrix(t(tau[gihj, ])), na.rm = TRUE) * Model$y2[g, ], na.rm = TRUE)
            }
            #denominator = sum(colSums(as.matrix(tau[gihj, ]), na.rm = TRUE) * Model$y2[g, ], na.rm = TRUE)
            denominator = 1 / (sd_xi ^ 2) + denominator
            x_seq = numerator / denominator
            x2_seq = x_seq ^ 2 + 1 / denominator
            return(list(x_seq = x_seq, x2_seq = x2_seq))
        }
    
    xx_loop <-
        function(i,
                 Model,
                 LFC,
                 x_seq_genes,
                 tau,
                 mean_xx,
                 sd_xx) {
            # Getting individual mean_xxij and sd_xxij values for each guide pair ij
            mean_xxij = mean_xx[i]
            sd_xxij = sd_xx[i]
            
            # check if both gihj and higj exist
            i_split = strsplit(x = i, split = Model$pattern_split, fixed = TRUE)[[1]]
            i_inv = paste0(rev(i_split), collapse = Model$pattern_split)
            if (!is.na(Model$xx[i_inv])){
                i = c(i,i_inv)
            }
            
            # Get genes gh and guides ij
            gi = Model$hashes_x$paired_guide[i, 1]
            hj = Model$hashes_x$paired_guide[i, 2]
            g = x_seq_genes[gi]
            h = x_seq_genes[hj]
            #gh = paste(sort(c(g, h)), collapse = Model$pattern_join)
            gh = vapply(seq_len(length(g)), function(x){
                return(paste(sort(c(g[x], h[x])), collapse = Model$pattern_join))
            }, character(1))
            
            numerator = Model$s[gh, ] * tau[i, ] * (LFC[i, ] - Model$x[gi] *
                                                        Model$y[g, ] - Model$x[hj] * Model$y[h, ])
            numerator = mean_xxij / (sd_xxij ^ 2) + sum(numerator, na.rm = TRUE)
            denominator = 1 / (sd_xxij ^ 2) + sum(Model$s2[gh, ] * tau[i, ], na.rm = TRUE)
            xx = numerator / denominator
            xx2 = xx ^ 2 + 1 / denominator
            return(list(xx = xx, xx2 = xx2))
        }
    if (verbose) {
        message("\tUpdating x:")
    }
    
    #### DEBUGGING: ###
    # x_res <- list()
    # for(x in names(Model$x)){
    #     x_res[[x]] <- try({
    #         x_loop(i = x,
    #            Model = Model,
    #            guide2gene = guide2gene,
    #            LFC = LFC,
    #            tau = tau,
    #            x_seq_genes = x_seq_genes,
    #            mean_x = mean_x,
    #            sd_x = sd_x)
    #     })
    #     if("try-error" %in% class(x_res[[x]])){
    #         print(x)
    #     }
    # }
    ###################
    group1 <- intersect(names(Model$x), setdiff(Model$hashes_x$paired_guide[,1], Model$hashes_x$paired_guide[,2]))
    group2 <- intersect(names(Model$x), setdiff(Model$hashes_x$paired_guide[,2], Model$hashes_x$paired_guide[,1]))
    group3 <- intersect(names(Model$x), intersect(Model$hashes_x$paired_guide[,1], Model$hashes_x$paired_guide[,2]))
    
    if(length(group1) > 0){
        if(verbose) {
            message("\t\tUpdating x for group 1...")
            x_group1 <- pbmcapply::pbmclapply(
                X = group1,
                FUN = x_loop,
                Model = Model,
                guide2gene = guide2gene,
                LFC = LFC,
                tau = tau,
                x_seq_genes = x_seq_genes,
                mean_x = mean_x,
                sd_x = sd_x,
                mc.cores = cores,
                ignore.interactive = verbose
            )
        }else{
            x_group1 <- parallel::mclapply(
                X = group1,
                FUN = x_loop,
                Model = Model,
                guide2gene = guide2gene,
                LFC = LFC,
                tau = tau,
                x_seq_genes = x_seq_genes,
                mean_x = mean_x,
                sd_x = sd_x,
                mc.cores = cores
            )
        }
        
        
        Model$x[group1] <- lapply(x_group1, magrittr::extract, "x_seq") %>%
            unlist(recursive = TRUE, use.names = FALSE)
        Model$x2[group1] <- lapply(x_group1, magrittr::extract, "x2_seq") %>%
            unlist(recursive = TRUE, use.names = FALSE)
    }

    if(length(group2) > 0){
        if(verbose){
            message("\t\tUpdating x for group 2...")
            x_group2 <- pbmcapply::pbmclapply(
                X = group2,
                FUN = x_loop,
                Model = Model,
                guide2gene = guide2gene,
                LFC = LFC,
                tau = tau,
                x_seq_genes = x_seq_genes,
                mean_x = mean_x,
                sd_x = sd_x,
                mc.cores = cores,
                ignore.interactive = verbose
            )
        }else{
            x_group2 <- parallel::mclapply(
                X = group2,
                FUN = x_loop,
                Model = Model,
                guide2gene = guide2gene,
                LFC = LFC,
                tau = tau,
                x_seq_genes = x_seq_genes,
                mean_x = mean_x,
                sd_x = sd_x,
                mc.cores = cores
            )
        }
        Model$x[group2] <- lapply(x_group2, magrittr::extract, "x_seq") %>%
            unlist(recursive = TRUE, use.names = FALSE)
        Model$x2[group2] <- lapply(x_group2, magrittr::extract, "x2_seq") %>%
            unlist(recursive = TRUE, use.names = FALSE)
    }
    
    if(length(group3) > 0){
        if(verbose){
            message("\t\tUpdating x for group 3...")
            pb <- pbmcapply::progressBar(max = length(group3))
        } 
        for(i in group3){
            i_res <- x_loop(Model = Model, 
                            i = i, 
                            guide2gene = guide2gene,
                            LFC = LFC,
                            tau = tau,
                            x_seq_genes = x_seq_genes,
                            mean_x = mean_x,
                            sd_x = sd_x)
            
            Model$x[i] <- magrittr::extract(i_res, "x_seq") %>%
                unlist(recursive = TRUE, use.names = FALSE)
            Model$x2[i] <- magrittr::extract(i_res, "x2_seq") %>%
                unlist(recursive = TRUE, use.names = FALSE)
            if(verbose)
                setTxtProgressBar(pb, value = getTxtProgressBar(pb)+1)
        }
        if(verbose)
            close(pb)
    }

    
    if (verbose) {
        message("\tUpdating xx:")
        xx_res <- pbmcapply::pbmclapply(
            names(Model$xx),
            FUN = xx_loop,
            Model = Model,
            LFC = LFC,
            x_seq_genes = x_seq_genes,
            tau = tau,
            mean_xx = mean_xx,
            sd_xx = sd_xx,
            mc.cores = cores,
            ignore.interactive = verbose
        )
    }else{
        xx_res <- parallel::mclapply(
            names(Model$xx),
            FUN = xx_loop,
            Model = Model,
            LFC = LFC,
            x_seq_genes = x_seq_genes,
            tau = tau,
            mean_xx = mean_xx,
            sd_xx = sd_xx,
            mc.cores = cores
            )
    }
    
    
    Model$xx[] <- lapply(xx_res, magrittr::extract, "xx") %>%
        unlist(use.names = FALSE)
    Model$xx2[] <- lapply(xx_res, magrittr::extract, "xx2") %>%
        unlist(use.names = FALSE)
    
    # output
    if (verbose) {
        tend = Sys.time()
        tdiff = difftime(tend, tstart)
        message("\tCompleted update of x.")
        message("\tTime to completion: ",
                round(tdiff, digits = 3),
                ' ',
                units(tdiff))
    }
    return(Model)
}
