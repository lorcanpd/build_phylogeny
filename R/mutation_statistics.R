

get_mutation_stats <- function(data){
    stats <- data %>%
    # First, annotate each mutation across all samples
    group_by(Muts) %>%
    mutate(
        # Count distinct samples with NV > 0 for each mutation
        sample_count = n_distinct(Sample[NV > 0]),
        # Flag as unique if only present in one sample
        is_unique = sample_count == 1
    ) %>%
    ungroup() %>%
    # Then, summarize per sample
    group_by(Sample) %>%
    summarise(
        total_mutations = n(),
        per_sample_mutations = sum(NV > 0),
        # Count unique mutations for the sample
        unique_mutations = sum(NV > 0 & is_unique),
        # Shared = Per sample - Unique
        shared_mutations = per_sample_mutations - unique_mutations
    )
    return(stats)
}