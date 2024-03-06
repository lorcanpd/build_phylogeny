
library(dplyr)
library(parallel)
library(doParallel)
library(VGAM)

create_germline_qval <- function(data, params, sex) {
    # This function calculates the germline q-value for each mutation in the
    # dataset. It uses a binomial test to evaluate the likelihood of observing
    # the given number of variant reads (NV) given the total number of reads
    # (NR) at that position across all samples. The q-values are adjusted for
    # multiple testing using the  Benjamini-Hochberg method. For samples with
    # copy number variations (CNVs), the q-value is set to 0 as these samples
    # are not suitable for this analysis.

    process_group <- function(subset) {
        NV_sum <- sum(subset$NV)
        NR_sum <- sum(subset$NR)
        has_CNV <- any(subset$has_CNV)

        germline_qval <- if (NR_sum == 0 || has_CNV) {
            # If the mutation is associated with CNV, assign q-value as 1.
            # This marks them as germline and they will be filtered out
            # later (and will be filtered if they are CNVs anyway).
            1
        } else {
            # If it's a non-CNV mutation in a female or male (autosomal), use a
            # binomial test with expected probability (p) of 0.5.
            if (sex == "female" ||
                (sex == "male" && !subset$is_XY_chromosomal[1])) {
                p.adjust(
                    binom.test(
                        x = NV_sum, n = NR_sum, p = 0.5, alt = 'less'
                    )$p.value,
                    method = "BH")
            # If it's a non-CNV mutation in a male (XY chromosomal), use a
            # binomial test with expected probability (p) of 0.95.
            } else if (sex == "male" && subset$is_XY_chromosomal[1]) {
                p.adjust(
                    binom.test(
                        x = NV_sum, n = NR_sum, p = 0.95, alt = 'less'
                    )$p.value,
                    method = "BH")
            } else {
                1  # Default case to handle any other scenarios.
            }
        }
        return(
            data.frame(
                Muts = unique(subset$Muts),
                germline_qval = germline_qval
            )
        )
    }

    # Split data by Muts
    data_split <- split(data, data$Muts)

    # Determine number of cores to use
    num_cores <- ifelse(is.null(params$ncores), 1, params$ncores)

    # Use mclapply for parallel processing if more than one core is specified
    if (num_cores > 1) {
        library(parallel)
        library(doParallel)
        cl <- makeCluster(params$ncores)
        registerDoParallel(cl)
        results <- mclapply(data_split, process_group, mc.cores = num_cores)
        stopCluster(cl)
    } else {
        # print(unique(data_split$Sample))
        # print(paste0("sum NV: ", sum(data_split$NV)))
        # print(paste0("sum NR: ", sum(data_split$NR)))
        results <- lapply(data_split, process_group)
    }

    # Combine results
    final_data <- do.call(rbind, results)

    # Join back to original data
    final_data <- left_join(data, final_data, by = "Muts")

    return(final_data)
}


estimate_beta_binomial_rho <- function(data, params) {
    # Function to estimate the overdispersion parameter for the beta-binomial
    # for each mutation. This is used to filter out mutations that are likely
    # to be artefacts. Computationally intensive, so will be parallelised if
    # more than one core is specified in the params file.

    # Define a grid of rho values to search. rho is bounded between 0 and 1, but
    # the beta-binomial model is undefined at 0 and 1, so we use a very small
    # value instead of 0, and 0.89 instead of 1.
    candidate_rho_values <- 10^seq(-6, -0.05, by=0.05)

    # Split data by mutation across all samples for processing
    data_split <- split(data, data$Muts)

    process_subset <- function(subset) {
        # subfunction to process a subset of the data, i.e. a single mutation
        # across all samples. If not parallelising, this function will be
        # called sequentially for each mutation. If parallelising, this
        # function will be called in parallel for each mutation.

        mean_VAF <- sum(subset$NV) / sum(subset$NR)
        log_likelihoods <- sapply(candidate_rho_values, function(rho) {
            sum(
                dbetabinom(
                    x = subset$NV, size = subset$NR, rho = rho, prob = mean_VAF,
                    log = TRUE
                )
            )
        })
        rho_estimate <- candidate_rho_values[which.max(log_likelihoods)]
        subset$rho_estimate <- rho_estimate
        return(subset)
    }

    # Determine number of cores to use
    num_cores <- ifelse(is.null(params$ncores), 1, params$ncores)

    # Use mclapply for parallel processing if more than one core is specified
    if (num_cores > 1) {
        library(parallel)
        library(doParallel)
        cl <- makeCluster(params$ncores)
        registerDoParallel(cl)
        results <- mclapply(data_split, process_subset, mc.cores = num_cores)
        stopCluster(cl)
    } else {
        results <- lapply(data_split, process_subset)
    }

    # Combine results and drop intermediate columns
    final_data <- bind_rows(results)

    return(final_data)
}

