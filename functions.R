
packages <- c("Rsamtools", "GenomicRanges", "dplyr")

for (package in packages) {
    if (!require(package, character.only = TRUE)) {
        library(package, character.only = TRUE)
    }
}

plot_spectrum <- function(
    bed, save, add_to_title="",
    genomeFile = "/nfs/cancer_ref01/Homo_sapiens/37/genome.fa"
){
    mutations <- as.data.frame(bed)
    colnames(mutations) <- c("chr", "pos", "ref", "mut")
    mutations$pos <- as.numeric(mutations$pos)
    mutations <- mutations[
        (mutations$ref %in% c("A", "C", "G", "T")) &
            (mutations$mut %in% c("A", "C", "G", "T")) &
            (mutations$chr %in% c(1:22, "X", "Y")),
    ]
    # Extract genomic ranges using given ranges determined by GRanges and
    # IRanges
    mutations$trinuc_ref <- as.vector(
        scanFa(
            genomeFile,
            GRanges(
                mutations$chr,
                IRanges(as.numeric(mutations$pos) - 1,
                        as.numeric(mutations$pos) + 1)
            )
        )
    )
    # 2. Annotating the mutation from the pyrimidine base
    ntcomp <- c(T="A", G="C", C="G", A="T")
    mutations <- mutations %>%
        mutate(
            sub = if_else(
                ref %in% c("A", "G"),
                paste(ntcomp[ref], ntcomp[mut], sep=">"),
                paste(ref, mut, sep=">")
            ),
            trinuc_ref_py = if_else(
                ref %in% c("A", "G"),
                paste(
                    ntcomp[rev(strsplit(trinuc_ref, split="")[[1]])],
                    collapse=""
                ),
                trinuc_ref
            )
        )
    # REPLACED WITH DPLYR MUTATE
    # mutations$sub <- paste(mutations$ref, mutations$mut, sep=">")
    # mutations$trinuc_ref_py <- mutations$trinuc_ref
    # for (j in 1:nrow(mutations)) {
    #     if (mutations$ref[j] %in% c("A","G")) { # Purine base
    #         mutations$sub[j] <- paste(ntcomp[mutations$ref[j]], ntcomp[mutations$mut[j]], sep=">")
    #         mutations$trinuc_ref_py[j] <- paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j], split="")[[1]])], collapse="")
    #     }
    # }


    # 3. Counting subs
    freqs <- table(
        paste(
            mutations$sub,
            paste(
                substr(mutations$trinuc_ref_py, 1, 1),
                substr(mutations$trinuc_ref_py, 3, 3),
                sep="-"
            ),
            sep=","
        )
    )
    sub_vec <- c("C>A","C>G","C>T","T>A","T>C","T>G")
    ctx_vec <- paste(
        rep(c("A","C","G","T"), each=4),
        rep(c("A","C","G","T"), times=4),
        sep="-"
    )
    full_vec <- paste(rep(sub_vec,each=16), rep(ctx_vec,times=6),sep=",")
    freqs_full <- freqs[full_vec]
    freqs_full[is.na(freqs_full)] <- 0
    names(freqs_full) <- full_vec

    xstr <- paste(
        substr(full_vec,5,5),
        substr(full_vec,1,1),
        substr(full_vec,7,7),
        sep=""
    )

    if(!is.null(save)) {
        pdf(save, width = 12, height = 4)
    }
    if(is.null(save)) {
        dev.new(width = 12, height = 4)
    }
    colvec <- rep(
        c("dodgerblue", "black", "red", "grey70", "olivedrab3", "plum2"),
        each=16
    )
    y <- freqs_full
    maxy <- max(y)
    h <- barplot(y, las=2, col=colvec, border=NA, ylim=c(0,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, ylab="Number mutations", main=paste0("Number of mutations: ",sum(freqs_full), add_to_title))
    for (j in 1:length(sub_vec)) {
        xpos <- h[c((j-1)*16+1,j*16)]
        rect(
            xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA,
            col=colvec[j*16]
        )
        text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
    }
    if(!is.null(save)) {
        dev.off()
    }
}


library(dplyr)

exact.binomial <- function(gender, NV, NR, cutoff=-5, qval_return=F) {
    # Function to filter out germline variants based on unmatched
    # variant calls of multiple samples from same individual (aggregate coverage
    # ideally >150 or so, but will work with less). NV is matrix of reads
    # supporting variants and NR the matrix with total depth (samples as
    # columns, mutations rows, with rownames as chr_pos_ref_alt or equivalent).
    # Returns a logical vector, TRUE if mutation is likely to be germline.

    # Convert NV and NR to data frames for dplyr processing
    NV_df <- as.data.frame(NV)
    NR_df <- as.data.frame(NR)

    # Combine NV and NR with their row names
    combined_df <- NV_df %>%
        mutate(chr_pos = rownames(NV_df)) %>%
        rowwise() %>%
        mutate(
            NV_sum = sum(c_across(-chr_pos)),
            NR_sum = sum(NR_df[chr_pos, ])
        )

    XY_chromosomal <- grepl("X|Y", combined_df$chr_pos)
    autosomal <- !XY_chromosomal

    # Perform binomial test
    combined_df <- combined_df %>%
        mutate(
            pval = case_when(
                # If condition is true, perform binomial test
                gender == "female" & NR_sum > 0 ~
                    binom.test(
                        x = NV_sum, n = NR_sum, p = 0.5, alt = 'less'
                    # reutn p-value
                    )$p.value,
                gender == "male" & autosomal & NR_sum > 0 ~
                    binom.test(
                        x = NV_sum, n = NR_sum, p = 0.5, alt = 'less'
                    )$p.value,
                gender == "male" & XY_chromosomal & NR_sum > 0 ~
                    binom.test(
                        x = NV_sum, n = NR_sum, p = 0.95, alt = 'less'
                    )$p.value,
                # If no conditions are met return 1 as default.
                TRUE ~ 1
            )
        )

    # Adjust p-values and determine germline status
    qval <- p.adjust(combined_df$pval, method = "BH")

    if (qval_return) {
        return(qval)
    } else {
        germline <- log10(qval) > cutoff
        return(germline)
    }
}

# REPLACED WITH DPLYR VERSION
# exact.binomial <- function(gender, NV, NR, cutoff=-5, qval_return=F) {
#     # Function to filter out germline variants based on unmatched
#     # variant calls of multiple samples from same individual (aggregate coverage
#     # ideally >150 or so, but will work with less). NV is matrix of reads
#     # supporting variants and NR the matrix with total depth (samples as
#     # columns, mutations rows, with rownames as chr_pos_ref_alt or equivalent).
#     # Returns a logical vector, TRUE if mutation is likely to be germline.
#
#     XY_chromosomal <- grepl("X|Y", rownames(NR))
#     autosomal <- !XY_chromosomal
#
#     if (gender == "female") {
#         NV_vec <- rowSums(NV)
#         NR_vec <- rowSums(NR)
#         pval <- rep(1, length(NV_vec))
#         for (n in 1:length(NV_vec)) {
#             if (NR_vec[n]>0){
#                 pval[n] <- binom.test(
#                     x=NV_vec[n],
#                     n=NR_vec[n],
#                     p=0.5,
#                     alt='less'
#                 )$p.value
#             }
#         }
#     }
#     # For male, split test in autosomal and XY chromosomal part
#     if (gender == "male") {
#         pval <- rep(1,nrow(NV))
#         NV_vec <- rowSums(NV)[autosomal]
#         NR_vec <- rowSums(NR)[autosomal]
#         pval_auto <- rep(1, sum(autosomal))
#         pval_XY <- rep(1, sum(XY_chromosomal))
#
#         for (n in 1:sum(autosomal)){
#             if (NR_vec[n]>0) {
#                 pval_auto[n] <- binom.test(
#                     x=NV_vec[n],
#                     n=NR_vec[n],
#                     p=0.5,
#                     alt='less'
#                 )$p.value
#             }
#         }
#         NV_vec <- rowSums(NV)[XY_chromosomal]
#         NR_vec <- rowSums(NR)[XY_chromosomal]
#         for (n in 1:sum(XY_chromosomal)) {
#             if (NR_vec[n] > 0) {
#                 pval_XY[n] <- binom.test(
#                     x=NV_vec[n],
#                     n=NR_vec[n],
#                     p=0.95,
#                     alt='less'
#                 )$p.value
#             }
#         }
#         pval[autosomal] <- pval_auto
#         pval[XY_chromosomal] <- pval_XY
#     }
#     qval <- p.adjust(pval, method="BH")
#     if (qval_return) {
#         return(qval)
#     } else {
#         germline <- log10(qval) > cutoff
#         return(germline)
#     }
# }


estimate_rho_gridml <- function(NV_vec, NR_vec) {
    # Function to estimate maximum likelihood value of rho for beta-binomial
    # model, using a grid search. NV_vec and NR_vec are vectors of variant
    # reads and total reads, respectively, for a single mutation across
    # multiple samples.

    # Define grid of rho values to search. rho is bounded between 0 and 1, but
    # the beta-binomial model is undefined at 0 and 1, so we use a very small
    # value instead of 0, and 0.89 instead of 1.
    rhovec <- 10^seq(-6, -0.05, by=0.05)
    # mu is the mean of the beta-binomial distribution, which is the same as
    # the mean of the binomial distribution, i.e. the mean of the VAFs.
    mu <- sum(NV_vec) / sum(NR_vec)
    # Calculate log-likelihood for each rho value
    ll <-  sapply(
        rhovec,
        function(rhoj) sum(dbetabinom(
            x=NV_vec, size=NR_vec, rho=rhoj, prob=mu, log=TRUE
        ))
    )
    # Return rho value with maximum log-likelihood
    return(rhovec[ll==max(ll)][1])
}

beta.binom.filter <- function(NR, NV) {
    # Function to apply beta-binomial filter for artefacts. Works best on sets
    # of clonal samples (ideally >10 or so). As before, takes NV and NR as
    # input. Optionally calculates pvalue of likelihood beta-binomial with
    # estimated rho fits better than binomial. This was supposed to protect
    # against low-depth variants, but use with caution. Returns logical vector
    # with good variants = TRUE

    # REPLACED WITH DPLYR VERSION
    # rho_est <- pval <- rep(NA,nrow(NR))
    # for (k in 1:nrow(NR)) {
    #     rho_est[k] <- estimate_rho_gridml(
    #         NV_vec = as.numeric(NV[k,]),
    #         NR_vec=as.numeric(NR[k,])
    #     )
    # }

    # Create a data frame from NR and NV for rowwise operation
    data <- data.frame(NR = I(NR), NV = I(NV))

    # Apply the estimate_rho_gridml function to each row
    data <- data %>%
        rowwise() %>%
        mutate(
            rho_est = estimate_rho_gridml(
                NV_vec = as.numeric(NV),
                NR_vec = as.numeric(NR)
            )
        ) %>%
        ungroup() # Ungroup to complete rowwise operations

    # Extract the rho_est column
    rho_est <- data$rho_est

    return(rho_est)
}


dbinomtrunc = function(x, size, prob, min_x=4) {
    # Function to calculate truncated binomial probability. x is the number of
    # successes, size is the number of trials, prob is the probability of
    # success, and minx is the minimum number of successes to consider.
    dbinom(x, size, prob) / pbinom(min_x-0.1, size, prob, lower.tail=FALSE)
}


# estep = function(x, size, p.vector, prop.vector, ncomp, mode) {
#     # Part of Expectation Maximisation algorithm for mixture models. x is the
#     # vector of variant reads, size is the vector of total reads, p.vector is
#     # the vector of probabilities for the individual components, prop.vector is
#     # the vector of proportions for the individual components, ncomp is the
#     # number of components, and mode is either "Truncated" or "Full" depending
#     # on whether the truncated binomial distribution is used or not.
#
#     p.mat_estep <- matrix(0, ncol=ncomp, nrow=length(x))
#     for (i in 1:ncomp) {
#         if (mode=="Truncated") {
#             p.mat_estep[, i] <- prop.vector[i] * dbinomtrunc(
#                 x, size, prob = p.vector[i]
#             )
#         }
#         if (mode=="Full") {
#             p.mat_estep[, i] <- prop.vector[i] * dbinom(
#                 x, size, prob = p.vector[i]
#             )
#         }
#     }
#
#     # Normalise the probabilities so that they sum to 1
#     norm <- rowSums(p.mat_estep)
#     p.mat_estep <- p.mat_estep / norm
#     # Calculate log-likelihood of the data given the model. Here we use the
#     # sum of the log of the probabilities, which is equivalent to the log of
#     # the product of the probabilities.
#     loglikelihood <- sum(log(norm))
#
#     # Here we assign each observation to the component with the highest
#     # probability. This is a crude way of assigning observations to clusters,
#     # and is only used for plotting purposes?
#     which_clust <- rep(1, length(x))
#     if (ncomp > 1) {
#         which_clust <- apply(p.mat_estep, 1, which.max)
#     }
#
#     return(
#         list(
#             "posterior" = p.mat_estep,
#             "LL" = loglikelihood,
#             "Which_cluster" = which_clust
#         )
#     )
# }

# Vectorised version of the estep function above should be faster using matrix
# operations instead of for loops.
# TODO test that this produces the same results as the original function and
#   is actually faster.
estep <- function(x, size, p.vector, prop.vector, ncomp, mode) {
    # Part of Expectation Maximisation algorithm for mixture models. x is the
    # vector of variant reads, size is the vector of total reads, p.vector is
    # the vector of probabilities for the individual components, prop.vector is
    # the vector of proportions for the individual components, ncomp is the
    # number of components, and mode is either "Truncated" or "Full" depending
    # on whether the truncated binomial distribution is used or not.

    # Create a matrix where each column is a repeat of 'x' and 'size'
    x_mat <- matrix(x, nrow = length(x), ncol = ncomp, byrow = TRUE)
    size_mat <- matrix(size, nrow = length(x), ncol = ncomp, byrow = TRUE)

    # Calculate the binomial probabilities for each component
    prob_mat <- if (mode == "Truncated") {
        dbinomtrunc(x_mat, size_mat, prob = p.vector)
    } else {
        dbinom(x_mat, size_mat, prob = p.vector)
    }

    # Multiply by proportions
    p.mat_estep <- prob_mat * prop.vector

    # Normalise the probabilities
    norm <- rowSums(p.mat_estep)
    p.mat_estep <- p.mat_estep / norm

    # Calculate log-likelihood
    loglikelihood <- sum(log(norm))

    # Assign each observation to the component with the highest probability
    which_clust <- apply(p.mat_estep, 1, which.max)

    return(
        list(
            "posterior" = p.mat_estep,
            "LL" = loglikelihood,
            "Which_cluster" = which_clust
        )
    )
}

mstep <- function(x, size, e.step) {
    # The maximisation step of the EM algorithm. x is the vector of variant
    # reads, size is the vector of total reads, and e.step is the output of
    # the estep function.

    # Calculate the mean of the posterior probabilities for each component of
    # the mixture model. The mean of these probabilities is interpreted as an
    # estimate of the proportion of the data that belongs to each component.
    prop.vector_temp <- colMeans(e.step$posterior)
    # Calcualte the weighted sum of the success probabilities for each component
    # of the mixture model. The weighted sum is weighted by the posterior
    # probabilities, and is interpreted as an estimate of the success
    # probability for each component. These are normalised by the sum of the
    # posterior probabilities for each component so that they sum to 1.
    p.vector_temp <- colSums(
        x / size * e.step$posterior) / colSums(e.step$posterior)

    list(
        "prop" = prop.vector_temp,
        "p" = p.vector_temp
    )
}

em.algo <- function(x, size, prop.vector_inits, p.vector_inits, maxit=5000,
                    tol=1e-6, nclust, binom_mode){
    ## prop.vector_inits =  initial values for the mixture proportions
    ## p.vector_inits =  initial values for the probabilities

    # Initiate EM
    flag <- 0
    e.step <- estep(x, size, p.vector = p.vector_inits, prop.vector = prop.vector_inits, ncomp=nclust, mode=binom_mode)
    m.step <- mstep(x, size, e.step)
    prop_cur <- m.step[["prop"]]
    p_cur <- m.step[["p"]]
    cur.LL <- e.step[["LL"]]
    LL.vector <- e.step[["LL"]]

    # Iterate between expectation and maximisation steps
    for (i in 2:maxit) {
        e.step <- estep(
            x, size, p.vector = p_cur, prop.vector = prop_cur, ncomp=nclust,
            mode=binom_mode
        )
        m.step <- mstep(x, size, e.step)
        prop_new <- m.step[["prop"]]
        p_new <- m.step[["p"]]

        prev_prop <- prop.vector_inits
        prev_p <- p.vector_inits

        LL.vector <- c(LL.vector, e.step[["LL"]])
        LL.diff <- abs((cur.LL - e.step[["LL"]]))
        which_clust <- e.step[["Which_cluster"]]

        # REPLACED WITH TRIPLE CHECK BELOW
        # Stop iteration if the difference between the current and new
        # log-likelihood is less than a tolerance level
        # if (LL.diff < tol) {
        #     flag <- 1
        #     break
        # }

        # Check for convergence in log-likelihood, mixture proportions, and
        # success probabilities. Checking the changes in the mixture proportions
        # and success probabilities is important because the log-likelihood
        # can appear converged, but the mixture proportions and
        # success probabilities can still oscillate - which is indicative of
        # a local maximum or plateau in the log-likelihood rather than
        # true convergence.
        if (sqrt(sum((prop_cur - prev_prop)^2)) < tol &&
            sqrt(sum((p_cur - prev_p)^2)) < tol &&
            LL_diff < tol) {
            flag <- 1
            break
        }
        prev_prop <- prop_cur
        prev_p <- p_cur
        # Otherwise continue iteration
        prop_cur <- prop_new
        p_cur <- p_new
        cur.LL <- e.step[["LL"]]
    }

    if (!flag) {
        warning("Did not converge\n")
    }

    BIC = log(length(x)) * nclust * 2 - 2 * cur.LL
    AIC = 4 * nclust - 2 * cur.LL
    list(
        "LL" = LL.vector,
        "prop" = prop_cur,
        "p" = p_cur,
        "BIC" = BIC,
        "AIC" = AIC,
        "n" = nclust,
        "Which_cluster" = which_clust
    )
}


binom_mix <- function(
    x, size, nrange=1:3, criterion="BIC", maxit=5000, tol=1e-6, mode="Full"
) {
    ## Perform the EM algorithm for different numbers of components
    ## Select best fit using the Bayesian Information Criterion (BIC)
    ## or the Akaike information criterion (AIC)
    i <- 1
    results <- list()
    BIC_vec <- c()
    AIC_vec <- c()

    # TODO: parallelise this loop using either the parallel package or
    #   foreach package. This is a good candidate for parallelisation because
    #   each iteration of the loop is independent of the others.
    for (n in nrange) {
        ## Initialise EM algorithm with values from kmeans clustering
        init <- kmeans(x / size, centers = n)
        prop_init <- init$size / length(x)
        p_init <- init$centers

        results[[i]] <- em.algo(
            x, size, prop.vector_inits = prop_init, p.vector_inits = p_init,
            nclust = n, maxit, tol, binom_mode = mode
        )
        BIC_vec <- c(BIC_vec, results[[i]]$BIC)
        AIC_vec <- c(AIC_vec, results[[i]]$AIC)
        i <- i + 1
    }
    # TODO come back to this later to see if the function could just return the
    #   results object and then the user can decide which criterion to use
    #   outside of the function.
    if (criterion=="BIC") {
        results[[which.min(BIC_vec)]]$BIC_vec <- BIC_vec
        return(results[[which.min(BIC_vec)]])
    }
    if (criterion=="AIC") {
        return(results[[which.min(AIC_vec)]])
    }
}

# Vectorised version of the binom_pval_matrix function below should be faster
# using matrix operations instead of for loops.
# TODO: test that this produces the same results as the original function and
#   is actually faster.
binom_pval_matrix <- function(NV, NR, gender, qval_return=FALSE) {
    # Replace zeros in NR with ones to avoid division by zero
    NR[NR == 0] <- 1

    # Define the success probability based on gender and chromosome
    p_success <- ifelse(
        gender == "male" & grepl("X|Y", rownames(NV)),
        0.95, 0.5
    )

    # Apply binomial test across the matrices
    pval_mat <- mapply(function(nv, nr, p) {
        if (nv > 0 && nr > 0) {
            binom.test(nv, nr, p, alternative = "less")$p.value
        } else {
            1  # Return 1 (or another appropriate value) for invalid or zero counts
        }
    }, NV, NR, p_success)

    # Reshape the output to a matrix with appropriate dimensions and names
    pval_mat <- matrix(pval_mat, nrow = nrow(NV), ncol = ncol(NV))
    rownames(pval_mat) <- rownames(NV)
    colnames(pval_mat) <- colnames(NV)

    # Adjust p-values if required
    if (qval_return) {
        qval_mat <- p.adjust(pval_mat, method = 'BH')
        return(qval_mat)
    } else {
        return(pval_mat)
    }
}

# REPLACED WITH VECTORISED VERSION ABOVE
# binom_pval_matrix <- function(NV, NR, gender, qval_return=FALSE) {
#     NR_nonzero <- NR
#     NR_nonzero[NR_nonzero==0] <- 1
#     pval_mat <- matrix(0, nrow = nrow(NV), ncol = ncol(NV))
#     rownames(pval_mat) <- rownames(NV)
#     colnames(pval_mat) <- colnames(NV)
#     if(gender == "male") {
#         for(i in 1:nrow(NV)) {
#             for (j in 1:ncol(NV)) {
#                 if (!grepl("X|Y",rownames(NV)[1])) {
#                     pval_mat[i, j] <- binom.test(
#                         NV[i, j], NR_nonzero[i, j], p = 0.5,
#                         alternative = "less"
#                     )$p.value
#                 }
#                 else {
#                     pval_mat[i, j] <- binom.test(
#                         NV[i, j], NR_nonzero[i, j], p = 0.95,
#                         alternative = "less"
#                     )$p.value
#                 }
#             }
#         }
#     } else if (gender == "female") {
#         for(i in 1:nrow(NV)) {
#             for (j in 1:ncol(NV)) {
#                 pval_mat[i, j] <- binom.test(
#                     NV[i, j], NR_nonzero[i, j], p = 0.5,
#                     alternative = "less"
#                 )$p.value
#             }
#         }
#     }
#     if (qval_return) {
#         qval_mat <- matrix(
#             p.adjust(as.vector(pval_mat), method='BH'),
#             ncol=ncol(pval_mat)
#         )
#         rownames(qval_mat) <- rownames(NV)
#         colnames(qval_mat) <- colnames(NV)
#         return(qval_mat)
#     } else {
#         return(pval_mat)
#     }
# }

# TODO: fix hardcoded path to output directory.
# TODO: File saving inside forloop is likely to make code very slow.
# TODO: Write error handling for cases without convergence etc.
# TODO: Try to vectorise and parallelise the for loop to make efficient.
# TODO: Refactor plotting to separate function for clarity and reusability.
# TODO: Add inline comments to explain each step of the function.
# TODO: Perhaps alter return of peak VAF? It may be more informative to return a
#   more comprehensive data structure that includes additional details about
#   the mixture model components for each sample.
apply_mix_model <- function(NV, NR, plot=TRUE, prop_cutoff=0.15) {
    peak_VAF <- rep(0,ncol(NV))
    names(peak_VAF) <- colnames(NV)
    autosomal <- !grepl("X|Y", rownames(NV))
    for(s in colnames(NV)) {
        muts_include <- NV[, s] > 3 & autosomal
        NV_vec <- NV[muts_include, s]
        NR_vec <- NR[muts_include, s]
        res <- binom_mix(NV_vec,NR_vec,mode="Truncated",nrange=1:3)
        saveRDS(res, paste0(output_dir, s, "_binom_mix.Rdata"))

        if (plot) {
            pdf(paste0(output_dir,s,"_binom_mix.pdf"))
            p <- hist(
                NV_vec/NR_vec, breaks=20, xlim=c(0, 1), col='gray', freq=FALSE,
                xlab="Variant Allele Frequency",
                main=paste0(s, ", (n=",length(NV_vec),")")
            )
            cols <- c("red","blue","green","magenta","cyan")

            y_coord <- max(p$density) - 0.5
            y_intv <- y_coord/5

            for (i in 1:res$n) {
                depth <- rpois(n=5000,lambda=median(NR_vec))
                sim_NV <- unlist(lapply(depth, rbinom, n=1, prob=res$p[i]))
                sim_VAF <- sim_NV/depth
                sim_VAF <- sim_VAF[sim_NV>3]
                dens <- density(sim_VAF)
                lines(
                    x=dens$x, y=res$prop[i]*dens$y, lwd=2, lty='dashed',
                    col=cols[i]
                )
                y_coord <- y_coord - y_intv / 2
                text(
                    y=y_coord, x=0.9,
                    label=paste0("p1: ", round(res$p[i], digits=2))
                )
                segments(
                    lwd=2, lty='dashed', col=cols[i], y0=y_coord + y_intv / 4,
                    x0=0.85, x1=0.95
                )
            }
            dev.off()
        }
        peak_VAF[s] <- max(res$p[res$prop > prop_cutoff])
    }
    return(peak_VAF)
}

add_ancestral_outgroup <- function(tree, outgroup_name="Ancestral") {
    # This function adds the ancestral tip at the end
    tmp <- tree$edge
    N <- length(tree$tip.label)
    newroot <- N + 2
    renamedroot <- N + 3
    ancestral_tip <- N + 1
    tmp <- ifelse(tmp > N, tmp + 2, tmp)

    tree$edge <- rbind(c(newroot, renamedroot), tmp, c(newroot,ancestral_tip))
    tree$edge.length <- c(0, tree$edge.length, 0)

    tree$tip.label <- c(tree$tip.label, outgroup_name)
    tree$Nnode <- tree$Nnode + 1
    mode(tree$Nnode) <- "integer"
    mode(tree$edge) <- "integer"
    return(tree)
}

