
library(dplyr)

# Functions for assigning flags to variants. These are used to optionally
# filter variants from the data.

determine_mutation_type <- function(data) {
    # Fucntion to determine the type of mutation (SNV, indel, or both).

    # I have recreated using the logic from the original script. However, I
    # think this could be simplified by removing the "both" category because
    # it is redundant. If the mutation is not an SNV or an indel, then it can't
    # be both?
    data <- data %>%
        mutate(
            Mutation_Type = case_when(
                nchar(Ref) > 1 | nchar(Alt) > 1 ~ "indel",
                nchar(Ref) == 1 & nchar(Alt) == 1 ~ "SNV",
                TRUE ~ "both"
            )
        )
    return(data)
}


create_XY_flag <- function(data) {
    # Function to determine whether a variant is autosomal or not.
    data$is_XY_chromosomal <- grepl("X|Y", data$Chr)
    return(data)
}


determine_sex <- function(data) {

    mean_autosomal_depth <- data %>%
        filter(!is_XY_chromosomal) %>%
        group_by(Sample) %>%
        summarise(per_sample_mean_depth = mean(NR, na.rm = TRUE)) %>%
        summarise(
            overall_mean_depth = mean(per_sample_mean_depth, na.rm = TRUE)
        )

    mean_XY_depth <- data %>%
        filter(is_XY_chromosomal) %>%
        group_by(Sample) %>%
        summarise(per_sample_mean_depth = mean(NR, na.rm = TRUE)) %>%
        summarise(
            overall_mean_depth = mean(per_sample_mean_depth, na.rm = TRUE)
        )

    if (mean_XY_depth > (0.8 * mean_autosomal_depth)) {
        return("female")
    } else {
        return("male")
    }
}


create_CNV_flag <- function(data, params) {
    # Create a flag to indicate whether a variant is in a sample with a CNV.
    data$has_CNV <- ifelse(
        data$Sample %in% params$samples_with_CNVs, TRUE, FALSE
    )
    return(data)
}


create_depth_flag <- function(data, params, sex) {
    data <- data %>%
        mutate(
            sufficient_depth = case_when(
                (!has_CNV &
                    sex == "male" &
                    NR > params$min_cov &
                    NR < params$max_cov &
                    !is_XY_chromosomal) ~ TRUE,
                (!has_CNV &
                    sex == "male" &
                    NR > params$min_cov / 2 &
                    NR < params$max_cov / 2 &
                    is_XY_chromosomal) ~ TRUE,
                (!has_CNV &
                    sex == "female" &
                    NR > params$min_cov &
                    NR < params$max_cov) ~ TRUE,
                TRUE ~ FALSE
            )
        )
    return(data)
}


# Sort data by column dplyir
sort_data <- function(data, column) {
    data <- data %>%
        arrange(column)
    return(data)
}

flag_shared_mutations <- function(data) {

    # Determine which mutations are shared
    shared_mutations <- data %>%
        filter(NV > 0) %>%
        group_by(Muts) %>%
        summarise(
            shared = n_distinct(Sample) >= 2,
            .groups = "drop"
        )

    # Join the shared flag back to the original data
    data <- data %>%
        left_join(shared_mutations, by = "Muts")

    # Ensure all rows have a defined shared value
    data$shared[is.na(data$shared)] <- FALSE

    return(data)
}

flag_conv_shared_mutations <- function(data, params) {

    # Determine which mutations are shared
    shared_mutations <- data %>%
        filter(NV > 0) %>%
        group_by(Muts) %>%
        summarise(
            conv_shared = n_distinct(Sample) >= params$min_variant_reads_shared,
            .groups = "drop"
        )

    # Join the shared flag back to the original data
    data <- data %>%
        left_join(shared_mutations, by = "Muts")

    # Ensure all rows have a defined shared value
    data$conv_shared[is.na(data$conv_shared)] <- FALSE

    return(data)
}
