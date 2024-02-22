
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


flag_shared_mutations <- function(data) {
    # Function to determine whether a mutation is shared between samples. This
    # is done by checking whether the mutation is present in at least two
    # samples. If it is, the mutation is flagged as shared. If not, it is
    # flagged as private.
    data <- data %>%
        group_by(Muts) %>%
        mutate(
            shared = ifelse(
                sum(NV > 0) > 1, TRUE, FALSE
            )
        )
    return(data)
}



