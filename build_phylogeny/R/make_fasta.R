
library(seqinr)
library(dplyr)
library(tidyr)


get_binom_pval <- function(data, sex) {

    process_sample <- function (sample_data, sex) {
        sample_data <- sample_data %>%
            rowwise() %>%
            mutate(
                binom_pval = case_when(
                    (
                        NV > 0 & NR > 0 &
                            (sex == "female" |
                                (sex == "male" &  !is_XY_chromosomal)
                            )
                    ) ~ binom.test(NV, NR, 0.5, alternative = "less")$p.value,
                    (
                        NV > 0 & NR > 0 &
                            sex == "male" & is_XY_chromosomal
                    ) ~ binom.test(NV, NR, 0.95, alternative = "less")$p.value,
                    TRUE ~ 1
                )
            ) %>%
            unnest(cols = c(binom_pval)) %>%
            mutate(
                binom_qval = p.adjust(binom_pval, method = "BH")
            ) %>%
            ungroup()

        return(sample_data)
    }

    # Split data by Sample
    data_split <- split(data, data$Sample)

    # Determine number of cores to use
    num_cores <- ifelse(is.null(params$ncores), 1, params$ncores)

    # Use mclapply for parallel processing if more than one core is specified
    if (num_cores > 1) {
        library(parallel)
        library(doParallel)
        cl <- makeCluster(params$ncores)
        registerDoParallel(cl)
        results <- mclapply(
            data_split, process_sample, sex, mc.cores = num_cores
        )
        stopCluster(cl)
    } else {
        results <- lapply(data_split, process_sample, sex)
    }
    # Combine results
    data <- do.call(rbind, results)

    return(data)
}

create_binary_genotype <- function(data, sex, params) {
    print("Generating almost-dichotomous genotype representations")
    if (params$only_snvs) {
        data <- data %>%
            filter(Mutation_Type == "SNV")
    }
    if (params$genotype_conv_prob) {
        data <- get_binom_pval(data, sex)
        if (!is.null(params$min_variant_reads_shared)) {
            data <- data %>%
                mutate(
                    has_sufficient_reads = ifelse(
                        !is.null(params$min_variant_reads_shared),
                        NV >= params$min_variant_reads_shared,
                        TRUE
                    ),
                    is_true_somatic = ifelse(
                        !is.null(params$min_pval_for_true_somatic_shared),
                        binom_pval > params$min_pval_for_true_somatic_shared,
                        TRUE
                    ),
                    has_sufficient_vaf = case_when(
                        (!is.null(params$min_vaf_shared) &
                            sex == "female") ~ VAF > params$min_vaf_shared,
                        (!is.null(params$min_vaf_shared) &
                            sex == "male" &
                            !is_XY_chromosomal) ~ VAF > params$min_vaf_shared,
                        (!is.null(params$min_vaf_shared) &
                            sex == "male" &
                            is_XY_chromosomal
                        ) ~ VAF > (params$min_vaf_shared * 2),
                        TRUE ~ TRUE
                    ),
                    binary_genotype = as.numeric(has_sufficient_reads &
                        is_true_somatic &
                        has_sufficient_vaf)
                ) %>%
                mutate(
                    # Here a value of 0.5 is assigned to mutations which are
                    # not confidently somatic, but have a p-value and number
                    # of reads which suggests that they may be somatic.
                    binary_genotype = case_when(
                        binary_genotype == 1 ~ 1,
                        (
                            binary_genotype == 0 & NV > 0 & binom_pval > 0.01
                        ) ~ 0.5,
                        (
                            binary_genotype == 0 & NV > 2 & binom_pval > 0.001
                        ) ~ 0.5,
                        (
                            binary_genotype == 0 & binom_pval > 0.05
                        ) ~ 0.5,
                        TRUE ~ 0
                    )
                ) %>%
                select(
                    Sample, Muts, binary_genotype, Mutation_Type
                )
        }
    } else {
        data <- data %>%
            mutate(
                binary_genotype = case_when(
                    (sex == "male" & !is_XY_chromosomal &
                        VAF < params$vaf_absent) ~ 0,
                    (sex == "male" & !is_XY_chromosomal &
                        VAF >= params$vaf_present) ~ 1,
                    (sex == "male" & is_XY_chromosomal &
                        VAF < (params$vaf_absent * 2)) ~ 0,
                    (sex == "male" & is_XY_chromosomal &
                        VAF >= (params$vaf_present * 2)) ~ 1,
                    (sex == "female" & VAF < params$vaf_absent) ~ 0,
                    (sex == "female" & VAF >= params$vaf_present) ~ 1,
                    TRUE ~ 0
                )
            ) %>%
            mutate(
                binary_genotype = case_when(
                    binary_genotype > 0 & binary_genotype < 1 ~ 0.5,
                    TRUE ~ binary_genotype
                )
            ) %>%
        select(Sample, Muts, binary_genotype, Mutation_Type)
    }

    return(data)
}

make_fasta <- function(binary_data, params) {
    # Function to generate a fasta file from a binary genotype matrix.

    dna_strings <- list()
    # Add inferred ancestor
    dna_strings[["Ancestral"]] <- paste(
        rep("A", length(unique(binary_data$Muts))),
        sep = "",
        collapse = ""
    )
    # Add samples
    for (sample in unique(binary_data$Sample)) {
        dna_strings[[sample]] <- paste(
            binary_data %>%
                filter(Sample == sample) %>%
                mutate(
                    bases = case_when(
                        binary_genotype == 0 ~ "A",
                        binary_genotype == 1 ~ "T",
                        binary_genotype == 0.5 ~ "?"
                    )
                ) %>%
                pull(bases),
            sep = "",
            collapse = ""
        )
    }

    # Write fasta file
    filename <- if (params$only_snvs) {
        paste0(params$output_dir, "/", params$donor_id, "_SNV_for_MPBoot.fa")
    } else {
        paste0(params$output_dir, "/", params$donor_id, "indel_for_MPBoot.fa")
    }
    write.fasta(
        dna_strings,
        names = names(dna_strings),
        file = filename
    )
}

