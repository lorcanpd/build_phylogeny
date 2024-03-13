

get_mutation_stats <- function(data) {

    if (!"shared" %in% colnames(data)) {
        # Flagging shared mutations
        data <- data %>%
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
            total_mutations = n(),
            mutations_with_supporting_read = sum(NV > 0),
            shared_mutations = sum(shared),
            unique_mutations = sum(NV > 0 & !shared)
        )

    return(stats)
}
