
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
    data <- data %>%
        filter(NV > 0)
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


get_cluster_flags <- function(sample_results, method = "BIC") {
    converged_results <- Filter(function(x) x$converged, sample_results)

    # Return an empty dataframe if no results converged
    if (length(converged_results) == 0) {
        return(data.frame(Which_cluster = integer(0), Row = integer(0)))
    }

    best_result_index <- if (method == "BIC") {
        which.min(sapply(converged_results, function(x) x$BIC))
    } else if (method == "AIC") {
        which.min(sapply(converged_results, function(x) x$AIC))
    } else {
        stop("Invalid method specified. Choose either 'BIC' or 'AIC'.")
    }

    best_result <- converged_results[[best_result_index]]

    if (!is.list(best_result) || !"Which_cluster" %in% names(best_result)) {
        stop("Best result does not contain 'Which_cluster' data.")
    }

    output <- best_result$Which_cluster %>%
        as.data.frame() %>%
        setNames("Which_cluster") %>%
        mutate(Row = row_number())

    return(output)
}


plot_all_mixture_models <- function(mixture_model_data, mixture_model_results, params) {
    cluster_flags <- lapply(mixture_model_results, get_cluster_flags)
    cluster_flags <- lapply(names(cluster_flags), function(sample_name) {
        cluster_flag <- cluster_flags[[sample_name]]
        cluster_flag$Sample <- sample_name
        return(cluster_flag)
    }) %>% bind_rows()

    mixture_model_data <- mixture_model_data %>%
        filter(NV > 0) %>%
        group_by(Sample) %>%
        mutate(Row = row_number()) %>%
        ungroup() %>%
        left_join(cluster_flags, by = c("Sample", "Row"))

    mixture_model_data$Which_cluster <- factor(mixture_model_data$Which_cluster, levels = 1:3)

    # Start PDF output
    pdf(paste0(params$output_dir, "/mixture_model_plots.pdf"), width = 12)

    # Loop over each sample and create plots
    for (sample in unique(mixture_model_data$Sample)) {
        sample_data <- mixture_model_data %>% filter(Sample == sample)

        # Create the stacked histogram plot
        p1 <- ggplot(sample_data, aes(x = VAF, fill = Which_cluster)) +
            geom_histogram(binwidth = 0.01) +
            geom_vline(aes(xintercept = peak_VAF), colour = "red", linetype = "dashed", size = 1) +
            theme_bw() +
            scale_fill_brewer(palette = "Set3") +
            theme(legend.position = "none") +
            ggtitle(paste("Stacked -", sample))

        # Create the overlapping histogram plot
        # p2 <- ggplot(sample_data, aes(x = VAF, fill = Which_cluster)) +
        #     # geom_histogram(binwidth = 0.01, position = "identity", alpha = 0.5) +
        #     geom_histogram(binwidth = 0.01, position = position_dodge(width = 0.01)) +
        #     geom_vline(aes(xintercept = peak_VAF), colour = "red", linetype = "dashed", size = 1) +
        #     theme_bw() +
        #     scale_fill_brewer(palette = "Dark2") +
        #     ggtitle(paste("Densities -", sample))
        p2 <- ggplot(sample_data, aes(x = VAF, fill = Which_cluster, color = Which_cluster)) +
            geom_density(alpha = 0.5, adjust = 1) +  # adjust parameter can be tweaked for smoothing
            geom_vline(aes(xintercept = peak_VAF), colour = "red", linetype = "dashed", size = 1) +
            theme_bw() +
            scale_fill_brewer(palette = "Set3") +
            scale_color_brewer(palette = "Set3") +  # Ensure outline colors match fill colors
            ggtitle(paste("Density -", sample)) +
            theme(legend.position = "right")  # Adjust legend position if needed


        # Arrange the two plots side by side for the current sample
        grid.arrange(p1, p2, ncol = 2)
    }

    # Close PDF output
    dev.off()
}
