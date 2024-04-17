

# Function that takes in our main mapping quality and base quality filtered
# mutation dataset and a dataset of the same mutations but with no quality
# filtering. It then calculates the ratio of supporting reads between the two
# and filters out mutations in the quality filtered dataset that have a ratio
# below a certain threshold.
# The intuition here is that if a greater number of reads support a mutation in
# the unfiltered dataset, then the read supporting the mutation are less likely
# to be real. This is because the ratio indicates that the region is one that
# is prone to poor mapping.

read_ratio_flags <- function (data, unfiltered_data, ratio_threshold = 0.75) {

    filtered <- data %>%
        select(Muts, Sample, NV) %>%
        dplyr::rename(NV_filtered = NV)

    unfiltered <- unfiltered_data %>%
        select(Muts, Sample, NV) %>%
        dplyr::rename(NV_unfiltered = NV)

    ratio_data <- filtered %>%
        left_join(unfiltered, by = c("Muts", "Sample")) %>%
        mutate(
            ratio = NV_filtered / NV_unfiltered,
#             filtered_by_read_ratio = ifelse(
#                 ratio < ratio_threshold, TRUE, FALSE
#             )
            filtered_by_read_ratio = case_when(
                is.na(ratio) ~ TRUE
                ratio < ratio_threshold ~ TRUE,
                ratio >= ratio_threshold ~ FALSE
                )
            )
        ) %>%
        select(Muts, Sample, filtered_by_read_ratio)

    # Join the read ratio flag back to the original data
    data <- data %>%
        left_join(ratio_data, by = c("Muts", "Sample"))

    return(data)
}


