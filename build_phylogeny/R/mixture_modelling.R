
dbinomtrunc <- function(x, size, prob, min_x=4) {
    # Function to calculate truncated binomial probability. x is the number of
    # successes, size is the number of trials, prob is the probability of
    # success, and minx is the minimum number of successes to consider.
    # dbinom(x, size, prob) / pbinom(min_x-0.1, size, prob, lower.tail=FALSE)
    denominator <- pbinom(min_x-0.1, size, prob, lower.tail=FALSE) + 1e-8
    dbinom(x, size, prob) / denominator
}

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
    prob_mat <- dbinomtrunc(x_mat, size_mat, prob = p.vector)
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


em.algo <- function(x, size, prop.vector_inits, p.vector_inits,
                    max_iterations=5000, tolerance=1e-6, nclust, binom_mode
){
    ## prop.vector_inits =  initial values for the mixture proportions
    ## p.vector_inits =  initial values for the probabilities

    # Initiate EM
    has_converged <- FALSE
    e.step <- estep(
        x, size, p.vector = p.vector_inits, prop.vector = prop.vector_inits,
        ncomp=nclust, mode=binom_mode
    )
    m.step <- mstep(x, size, e.step)
    prop_cur <- m.step[["prop"]]
    p_cur <- m.step[["p"]]
    cur.LL <- e.step[["LL"]]
    LL.vector <- e.step[["LL"]]

    # Iterate between expectation and maximisation steps
    for (i in 2:max_iterations) {
        e.step <- estep(
            x, size, p.vector = p_cur, prop.vector = prop_cur, ncomp=nclust,
            mode=binom_mode
        )
        m.step <- mstep(x, size, e.step)
        prop_new <- m.step[["prop"]]
        p_new <- m.step[["p"]]

        # prev_prop <- prop.vector_inits
        # prev_p <- p.vector_inits

        LL.vector <- c(LL.vector, e.step[["LL"]])
        LL_diff <- abs((cur.LL - e.step[["LL"]]))
        which_clust <- e.step[["Which_cluster"]]

        # Check for convergence in log-likelihood, mixture proportions, and
        # success probabilities. Checking the changes in the mixture proportions
        # and success probabilities is important because the log-likelihood
        # can appear converged, but the mixture proportions and success
        # probabilities can still oscillate - which is indicative of a local
        # maximum or plateau in the log-likelihood rather than true convergence.

        if (LL_diff < tolerance) {
            has_converged <- TRUE
            break
        }
        # prev_prop <- prop_cur
        # prev_p <- p_cur
        # Otherwise continue iteration
        prop_cur <- prop_new
        p_cur <- p_new
        cur.LL <- e.step[["LL"]]
    }

    if (!has_converged) {
        warning("Did not converge\n")
    }

    BIC = log(length(x)) * nclust * 2 - 2 * cur.LL
    AIC = 4 * nclust - 2 * cur.LL
    list(
        # "LL" = LL.vector,
        "prop" = prop_cur,
        "p" = p_cur,
        "BIC" = BIC,
        "AIC" = AIC,
        "n" = nclust,
        "Which_cluster" = which_clust,
        "converged" = has_converged
    )
}


binomial_mixture <- function(
    x, size, nrange=1:3, criterion="BIC", max_iterations=5000,
    tolerance=1e-6, mode="Full"
) {
    ## Perform the EM algorithm for different numbers of components
    ## Select best fit using the Bayesian Information Criterion (BIC)
    ## or the Akaike information criterion (AIC)
    i <- 1
    results <- list()
    # BIC_vec <- c()
    # AIC_vec <- c()

    num_cores <- ifelse(is.null(params$ncores), 1, params$ncores)

    if (num_cores > 1) {
        library(parallel)
        library(doParallel)
        cl <- makeCluster(params$ncores)
        registerDoParallel(cl)

        results <- mclapply(
            nrange, function(n) {
                ## Initialise EM algorithm with values from kmeans clustering
                init <- kmeans(x / size, centers = n)
                prop_init <- init$size / length(x)
                p_init <- init$centers
                em.algo(
                    x, size, prop.vector_inits = prop_init,
                    p.vector_inits = p_init, nclust = n, max_iterations,
                    tolerance, binom_mode = mode
                )
        }, mc.cores = num_cores)

        # BIC_vec <- sapply(results, function(r) r$BIC)
        # AIC_vec <- sapply(results, function(r) r$AIC)

        stopCluster(cl)
    } else {
        for (n in nrange) {
            init <- kmeans(x / size, centers = n)
            prop_init <- init$size / length(x)
            p_init <- init$centers
            result <- em.algo(
                x, size, prop.vector_inits = prop_init, p.vector_inits = p_init,
                nclust = n, max_iterations, tolerance, binom_mode = mode
            )
            results[[n]] <- result
            # BIC_vec[n] <- result$BIC
            # AIC_vec[n] <- result$AIC
        }
    }

    return(results)
}

fit_binom_mix_model <- function(
    data, params, tolerance=1e-6, max_iteratons=5000, prop_cutoff=0.15
) {
    # Fit a binomial mixture model to the data.
    data_split <- split(data, data$Sample)

    process_sample <- function (sample_data) {
        # Helper function to process each sample in parallel and return the
        # peak VAF for each sample.
        binom_results <- binomial_mixture(
            sample_data$NV, sample_data$NR, tolerance = tolerance,
            max_iterations = max_iteratons, mode="Truncated", nrange=1:3
        )
    }

    num_cores <- ifelse(is.null(params$ncores), 1, params$ncores)

    # Use mclapply for parallel processing if more than one core is specified
    if (num_cores > 1) {
        library(parallel)
        library(doParallel)
        cl <- makeCluster(params$ncores)
        registerDoParallel(cl)
        results <- mclapply(data_split, process_sample, mc.cores = num_cores)
        stopCluster(cl)
    } else {
        results <- lapply(data_split, process_sample)
    }

    return(results)
}

process_mixture_model_results <- function(sample_results, method = "BIC") {
    # Filter results for those where convergence is TRUE
    converged_results <- Filter(function(x) x$converged, sample_results)

    # If there are no converged results, return NA
    if (length(converged_results) == 0) {
        return(data.frame(peak_VAF = NA, n_components = NA, BIC = NA))
    }

    if (method == "BIC") {
            # Select the result with the minimum BIC
    best_result <- converged_results[[
        which.min(sapply(converged_results, function(x) x$BIC))
    ]]
    } else if (method == "AIC") {
        # Select the result with the minimum number of components
        best_result <- converged_results[[
            which.min(sapply(converged_results, function(x) x$AIC))
        ]]
    }

    # Extract the peak VAF from the best result
    # Adjust the prop threshold as needed
    peak_VAF <- max(best_result$p[best_result$prop > 0.15])
    n_components <- best_result$n

    if (method == "BIC") {
        BIC <- best_result$BIC
        results <- data.frame(
            peak_VAF = peak_VAF, n_components = n_components, BIC = BIC
        )
    } else if (method == "AIC") {
        AIC <- best_result$AIC
        results <- data.frame(
            peak_VAF = peak_VAF, n_components = n_components, AIC = AIC
        )
    }

    return(results)
}

get_cluster_flags <-function(sample_results, method = "BIC"){
    converged_results <- Filter(function(x) x$converged, sample_results)

    # If there are no converged results, return NA
    # if (length(converged_results) == 0) {
    #     return(data.frame(peak_VAF = NA, n_components = NA, BIC = NA))
    # }

    if (method == "BIC") {
            # Select the result with the minimum BIC
    best_result <- converged_results[[
        which.min(sapply(converged_results, function(x) x$BIC))
    ]]
    } else if (method == "AIC") {
        # Select the result with the minimum number of components
        best_result <- converged_results[[
            which.min(sapply(converged_results, function(x) x$AIC))
        ]]
    }

    # Convert which cluster to a long format table with sample names row names
    # and which cluster as a column
    output <- best_result$Which_cluster %>%
        as.data.frame() %>%
        setNames("Which_cluster") %>%
        mutate(Row = row_number())


    return(output)
}

get_peak_VAFs <- function (mixture_model_data, mixture_model_results) {
    # Get the peak VAF for each sample from the mixture model results
    peak_VAFs <- lapply(mixture_model_results, process_mixture_model_results)

    # Convert to a data frame
    peak_VAFs_df <- do.call(rbind, peak_VAFs)
    # Make rownames a column
    peak_VAFs_df <- peak_VAFs_df %>%
        rownames_to_column("Sample")

    # Add the peak VAF to the mixture model data
    mixture_model_data <- mixture_model_data %>%
        left_join(
            peak_VAFs_df,
            by = c("Sample")
        )

    return(mixture_model_data)
}

plot_all_mixture_models <- function(
    mixture_model_data, mixture_model_results, params
) {
    # Plot a histogram of the VAF data from mixture_model_data and colour the plot
    # using the mixture model results which_cluster column
    # First assign the which cluster column to the mixture_model_data
    cluster_flags <- lapply(mixture_model_results, get_cluster_flags)


    cluster_flags <- lapply(
        names(cluster_flags), function(sample_name) {
            # Add the sample name as a new column
            cluster_flags[[sample_name]] %>%
            mutate(Sample = sample_name)
        }) %>%
        bind_rows()

    mixture_model_data <- mixture_model_data %>%
        group_by(Sample) %>%
        mutate(Row = row_number()) %>%
        ungroup() %>%
        left_join(
            cluster_flags,
            by = c("Sample", "Row")
        )

    # For each sample, plot a histogram of the VAF data from mixture_model_data and
    # colour the plot using the mixture model results Which_cluster column, add a
    # vertical line at the peak VAF.
    library(ggplot2)
    pdf(paste0(params$output_dir, "mixture_model_plots.pdf"))
    mixture_model_data %>%
        ggplot(aes(x = VAF)) +
        # colour by which cluster
        geom_histogram(
            data = . %>% filter(!is.na(Which_cluster)),
            aes(fill = Which_cluster),
            binwidth = 0.01,
            colour = "black"
        ) +
        geom_vline(
            aes(xintercept = peak_VAF),
            colour = "red",
            linetype = "dashed",
            size = 1
        ) +
        facet_wrap(~ Sample) +
        theme_bw()
    dev.off()
}

# OLD CODE GENERATES DATA TO CHECK THAT MODEL IS WORKING CORRECTLY, WHICH SEEMS
# WRONG. SHOULD BE USING REAL DATA TO CHECK THAT THE MODEL IS WORKING CORRECTLY.
# plot_binom_mix_model <- function(binom_mix_results, params) {
#
#     plot_binom_mix_model_sample <- function(sample_data) {
#         p <- hist(
#             sample_data$VAF, breaks = 20, plot = FALSE, xlab = "VAF",
#             xlim = c(0, 1),
#             col="grey", freq = FALSE, xlabel = "Variant Allele Frequency",
#             main = paste0(
#                 "Sample:", unique(sample_data$Sample), " (n=", nrow(data), ")"
#             )
#         )
#         cols <- c("red", "blue", "green", "magenta", "cyan")
#
#         y_coord <- max(p$density) - 0.5
#         y_intv <- y_coord / 5
#
#         for (i in 1:binom_mix_results$n) {
#             depth <- rpois(1, median(sample_data$NR))
#             simulated_NV <- rbinom(length(depth), depth, binom_mix_results$p[i])
#             simulated_VAF <- simulated_NV[simulated_NV > 0] / depth
#             dens <- density(simulated_VAF)
#             lines(
#                 dens$x, binom_mix_results$prop[i] * dens$y,  col = cols[i],
#                 lwd = 2, lty = "dashed"
#             )
#             text(
#                 y=y_coord, x=0.9,
#                 labels=paste0("p1: ",round(res$p[i],digits=2))
#             )
#             segments(
#                 lwd=2, lty="dashed", col=cols[i], y0=y_coord + y_intv / 4,
#                 x0=0.85, x1=0.95
#             )
#         }
#         dev.off()
#     }
#     num_cores <- ifelse(is.null(params$ncores), 1, params$ncores)
# }




