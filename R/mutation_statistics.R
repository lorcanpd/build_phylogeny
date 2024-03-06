

get_mutation_stats <- function(data) {

    if (!"shared" %in% colnames(data)) {
        # Flagging shared mutations
        data <- data %>%
            filter(NV > 0) %>%
            group_by(Muts) %>%
            mutate(
                # Identify shared mutations as those present in more than one sample
                shared = n_distinct(Sample) > 1
            ) %>%
            ungroup()
    }


    # Calculate stats per sample with the newly flagged shared status
    stats <- data %>%
        group_by(Sample) %>%
        summarise(
            total_mutations = n_distinct(Muts),
            mutations_with_supporting_read = sum(NV > 0),
            shared_mutations = sum(shared),
            unique_mutations = sum(NV > 0 & !shared),
            .groups = "drop"
        )

    return(stats)
}
