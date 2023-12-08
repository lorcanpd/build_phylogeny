
library(extraDistr)

set.seed(123) # Setting a seed for reproducibility

create_unique_mutations <- function(num_mutations, sex, y_ratio) {
    # This function creates a vector of mutations, where each mutation is
    # represented as a string in the format "chr_pos_ref_alt". The number of
    # mutations on the y chromosome is determined by the y_ratio parameter.

    mutations <- vector("character", num_mutations)
    chromosomes <-
        rep(paste0("chr", 1:22), times = 10)

    if (sex == "male") {
        # Adjust the frequency of chrY based on the y_ratio
        chromosomes <- c(
            chromosomes,
            rep("chrY", times = round(10 * y_ratio)),
            rep("chrX", times = 5)
        )
    } else {
        rep("chrX", times = 10)
    }

    for (i in 1:num_mutations) {
        chromosome <- sample(chromosomes, 1)
        position <- sample(100000000:999999999, 1)
        ref <- sample(c("A", "T", "C", "G"), 1)
        # Ensure alt is different from ref
        alt <- sample(setdiff(c("A", "T", "C", "G"), ref), 1)
        mutations[i] <- paste0(chromosome, "_", position, "_", ref, "_", alt)
    }
    return(mutations)
}

simulate_artefacts <- function(n, size) {
    # n: number of artefacts to simulate
    # size: total read counts (NR values)

    alpha <- rlnorm(n, meanlog = log(0.4), sdlog = sqrt(0.2))
    beta <- rlnorm(n, meanlog = log(0.4), sdlog = sqrt(0.2))
    # Simulate artefact variant read counts (NV) using beta-binomial distribution
    rbbinom(n, size, alpha, beta)
}

generate_nv_nr_data <- function(
    sex, num_mutations = 2560, num_samples = 320, avg_coverage = 30,
    y_ratio = 0.2, num_subclones = 2, germline_proportion = 0.5,
    artefact_proportion = 0.3
) {
    # Generate sample and mutation names
    samples <- paste0("Sample_", seq(1, num_samples))
    mutations <- create_unique_mutations(num_mutations, sex, y_ratio)

    # Step 1: Determine Germline Mutations, artefact mutations, and autosomal
    germline_indices <- sample(num_mutations, num_mutations * germline_proportion)
    is_germline <- rep(FALSE, num_mutations)
    is_germline[germline_indices] <- TRUE

    num_artefacts <- round(length(is_germline == FALSE) * artefact_proportion)
    artefact_indices <- sample(num_mutations, num_artefacts)

    is_autosomal <- !grepl("chrX|chrY", mutations)

    # Step 2: Assign Subclones to Samples
    sample_subclone_proportions <- matrix(
        runif(num_samples * num_subclones, min = 0, max = 1),
        nrow = num_samples, ncol = num_subclones
    )

    # Ensure that the sum of subclone proportions for each sample is 0.5
    sample_subclone_proportions <- t(apply(
        sample_subclone_proportions, 1, function(x) {
            x / sum(x) * 0.5
        }
    ))

    # Step 3: Assign Mutations to Subclones
    mutation_subclone_assignment <- sample(
        num_subclones, num_mutations, replace = TRUE
    )

    # Generate NR matrix
    NR <- matrix(
        rpois(num_mutations * num_samples, lambda = avg_coverage),
        nrow = num_mutations, ncol = num_samples,
        dimnames = list(mutations, samples)
    )
    NR[grepl("chrY", mutations)] <- round(NR[grepl("chrY", mutations)] * 0.4)
    NR[grepl("chrX", mutations)] <- round(NR[grepl("chrX", mutations)] * 0.5)

    # Generate NV matrix
    NV <- matrix(0, nrow = num_mutations, ncol = num_samples)
    colnames(NR) <- colnames(NV) <- samples
    rownames(NR) <- rownames(NV) <- mutations

    # Step 4: Sample Read Counts Based on Subclone Proportions
    NV <- apply(NR, 2, function(nr_vector, sample_index) {
        nv_vector <- numeric(num_mutations)
        for (i in 1:num_mutations) {
            if (is_germline[i]) {
                mean_proportion <- ifelse(
                    grepl("chrX|chrY", mutations[i]), 0.95, 0.5
                )
                nv_vector[i] <- rbinom(1, nr_vector[i], mean_proportion)
            } else if (i %in% artefact_indices) {
                # Artefact mutation simulation
                nv_vector[i] <- simulate_artefacts(1, nr_vector[i])
            } else {
                subclone_id <- mutation_subclone_assignment[i]
                if (sex == "male" & !is_autosomal[i]) {
                    subclone_proportion <- sample_subclone_proportions[
                        sample_index, subclone_id
                    ] * 2
                } else {
                    subclone_proportion <- sample_subclone_proportions[
                        sample_index, subclone_id
                    ]
                }
                nv_vector[i] <- rbinom(1, nr_vector[i], subclone_proportion)
            }
        }
        return(nv_vector)
    }, sample_index = 1:num_samples)
    rownames(NV) <- mutations
    return(list(NV = NV, NR = NR))
}



# # TODO use this function to generate data for unit tests
# # EXAMPLE OF USAGE.
# # Generate NV and NR data
# data <- generate_nv_nr_data(
#     sex = "male", num_mutations = 2560, num_samples = 8, avg_coverage = 30,
#     y_ratio = 0.2
# )

data <- generate_nv_nr_data(
    sex = "male", num_mutations = 5120, num_samples = 4, avg_coverage = 30,
    y_ratio = 0.2, num_subclones = 2, germline_proportion = 0.70,
    artefact_proportion = 0.2
)

NV <- data$NV
NR <- data$NR

# Display the first few rows of NV and NR
head(NV)
head(NR)

# Save NV and NR matrices to txt file
write.table(NV, file = "unit_tests/NV_dummy.txt", sep = "\t", quote = FALSE)
write.table(NR, file = "unit_tests/NR_dummy.txt", sep = "\t", quote = FALSE)


