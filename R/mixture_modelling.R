
mixture_modelling <- function(
    data, params, tolerance = 1e-6, max_iter = 5000, method = "BIC") {

    # Split data by Sample
    data_split <- split(data, data$Sample)

    process_sample <- function(sample_data) {
        # Get NV and NR columns plus muts as rowname in matrix
        x <- as.matrix(sample_data[, c("NV", "NR")])
        rownames(x) <- sample_data$Muts

        models <- list()
        aic_scores <- numeric()
        bic_scores <- numeric()

        for (n_clusters in 1:3) {
            model <- flexmix(
                x ~ 1, k = n_clusters,
                model = FLXMRglm(family = "binomial"),
                control = list(iter.max = max_iter, tolerance = tolerance)
            )
            models[[n_clusters]] <- model
            aic_scores[n_clusters] <- AIC(model)
            bic_scores[n_clusters] <- BIC(model)
        }

        # Determine the best model based on the lowest AIC score
        if (method == "BIC") {
            best_model_index <- which.min(bic_scores)
        } else if (method == "AIC") {
            best_model_index <- which.min(aic_scores)
        } else {
            stop("Invalid method specified. Choose either 'BIC' or 'AIC'.")
        }
        best_model <- models[[best_model_index]]

        # Cluster proportions
        cluster_proportions <- best_model@prior

        # Parameter estimates for each component
        probs <- c()
        for (n in 1:length(best_model@components)) {
            # Extract the log odds
            log_odds <- best_model@components[[n]][[1]]@parameters$coef
            # Convert to probabilities
            probs <- c(probs, 1 / (1 + exp(-log_odds)))
        }
        # Peak VAF
        peak_VAF <- probs[which.max(cluster_proportions)]

        # Join cluster assignments, peak VAF, proportion of largest cluster to
        # their respective Muts id.
        clusters <- clusters(best_model)
        clusters <- data.frame(clusters)
        rownames(clusters) <- rownames(x)
        # Convert rownames to a column
        clusters <- clusters %>%
            rownames_to_column("Muts") %>%
            setNames(c("Muts", "cluster"))


        # join the cluster assignments to the original data
        sample_data <- sample_data %>%
            left_join(
                clusters,
                by = c("Muts")
            ) %>%
            mutate(
                peak_VAF = peak_VAF
            )

        return(sample_data)
    }

    # Determine number of cores to use
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

    data <- bind_rows(results)

    return(data)
}


plot_mixture_models <- function(data, params) {
    pdf(paste0(params$output_dir, "/mixture_model_plots.pdf"), width = 11.7, height = 8.3)

    cluster_palette <- c("1" = "#8DD3C7", "2" = "#FFFFB3", "3" = "#BEBADA")
    line_palette <- c(
        "Subclone VAF threshold" = "#FB8072",
        "Main subclone peak VAF" = "#80B1D3"
    )

    for (sample in unique(data$Sample)) {
        sample_data <- data %>% filter(Sample == sample)
        peak_VAF <- unique(sample_data$peak_VAF)

        # Plot clusters with peak VAF line
        p1 <- ggplot(sample_data, aes(x = VAF, fill = as.factor(cluster))) +
            geom_histogram(bins = 30, position = "stack", show.legend = FALSE) +
            geom_vline(aes(xintercept = peak_VAF, color = "Main subclone peak VAF"), linetype = "dashed", size = 1.5) +
            geom_vline(aes(xintercept = params$vaf_threshold_mixmodel, color = "Subclone VAF threshold"), linetype = "dotted", size = 1.5) +
            scale_fill_manual(values = cluster_palette) +
            scale_color_manual(values = line_palette, name = "Thresholds") +
            # guides(fill = guide_legend(title = "Clusters", order = 1)) +
            theme_bw() +
            theme(legend.position = "bottom") +
            ggtitle(paste("Stacked -", sample)) +
            coord_cartesian(xlim = c(0, 1))

        # Density Plot
        p2 <- ggplot(sample_data, aes(x = VAF, fill = as.factor(cluster))) +
            geom_density(alpha = 0.7, adjust = 1, color = "black") +
            geom_vline(aes(xintercept = peak_VAF, color = "Main subclone peak VAF"), linetype = "dashed", size = 1.5, show.legend=FALSE) +
            geom_vline(aes(xintercept = params$vaf_threshold_mixmodel, color = "Subclone VAF threshold"), linetype = "dotted", size = 1.5, show.legend=FALSE) +
            scale_fill_manual(values = cluster_palette, name = "Subclone clusters") +
            scale_color_manual(values = line_palette, name = "") +
            theme_bw() +
            theme(legend.position = "bottom") +
            ggtitle(paste("Density -", sample)) +
            coord_cartesian(xlim = c(0, 1))

        # Arrange the two plots side by side for the current sample
        print(grid.arrange(p1, p2, ncol = 2))
    }
    dev.off()
}

