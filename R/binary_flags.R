

flag_binary_germline <- function(data, presence_threshold) {
    # If there are at least N supporting reads for a mutation across all the
    # samples we flag it with binary_germline == TRUE.
    data <- data %>%
        group_by(Muts) %>%
        mutate(binary_germline = all(NV >= presence_threshold)) %>%
        ungroup()

    return(data)
}


flag_binary_shared <- function(data, presence_threshold) {
    data <- data %>%
        group_by(Muts) %>%
        mutate(binary_shared = sum(NV >= presence_threshold) > 1) %>%
        ungroup()

    return(data)
}


flag_binary_unique <- function(data, presence_threshold) {
    data <- data %>%
        group_by(Muts) %>%
        mutate(binary_unique = sum(NV >= presence_threshold) == 1) %>%
        ungroup()

    return(data)
}
