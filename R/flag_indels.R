

flag_close_to_indel <- function (mutations, params) {

    process_subset <- function (mutations_subset, indel_subset) {
        # Adjusting indel_subset to include the window
        indel_subset <- indel_subset %>%
            mutate(
                start = Pos - params$indel_window,
                end = Pos + nchar(Ref) + params$indel_window
            ) %>%
            # As splitting by chromosome we don't need to keep the chromosome,
            # we just need the indel regions.
            select(start, end)


        # Check for overlap using a non-equi join
        result <- mutations_subset %>%
            rowwise() %>%
            mutate(
                close_to_indel = any(
                    Pos >= indel_subset$start &
                    Pos <= indel_subset$end
                )
            ) %>%
            select(Muts, close_to_indel)

        return(result)
    }

    indels <- fread(
        params$indel_file, header = FALSE, sep = "\t", data.table = TRUE
    )
    colnames(indels) <- c("Chr", "Pos", "Ref", "Alt", "Qual")
    mutation_locations <- mutations %>%
        select(Muts, Chr, Pos) %>%
        group_by(Muts, Chr, Pos) %>%
        distinct() %>%
        ungroup()

    # Split data and indel file by chromosome.
    mutations_by_chr <- split(mutation_locations, mutation_locations$Chr)
    indel_by_chr <- split(indels, indels$Chr)
    # Determine number of cores to use
    num_cores <- ifelse(is.null(params$ncores), 1, params$ncores)

    # Use mclapply for parallel processing if more than one core is specified
    if (num_cores > 1) {
        library(parallel)
        library(doParallel)
        cl <- makeCluster(params$ncores)
        registerDoParallel(cl)
        results <- mclapply(names(mutations_by_chr), function(chr) {
            process_subset(mutations_by_chr[[chr]], indel_by_chr[[chr]])
        }, mc.cores = num_cores)
        stopCluster(cl)
    } else {
        results <- lapply(names(mutations_by_chr), function(chr) {
            process_subset(mutations_by_chr[[chr]], indel_by_chr[[chr]])
        })
    }

    # Combine results and join flag back to original data
    results <- bind_rows(results)

    mutations <- mutations %>%
        left_join(results, by = c("Muts"))

    return(mutations)
}
