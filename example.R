
# This script showcases how the modularised pipeline can be used to build a
# phylogenetic tree, with optional iterative filtering and plotting of
# mutational signatures. The pipeline is designed to be flexible and modular,
# allowing for easy modification and extension. Uses are expected to write a
# custom pipeline which uses the functionality demonstrated here to suit their
# specific needs.

library(buildphylogeny)

print("Reading data")
params <- read_json_params("example_params.json")
data <- process_data(params)
data <- annotate_mutations(data, genomeFile = params$genome_file)

print("Generating filtering flags")
data <- create_XY_flag(data)
sex <- determine_sex(data)
data <- create_CNV_flag(data, params)
data <- determine_mutation_type(data)
data <- create_depth_flag(data, params, sex)

print("Estimating germline q-values")
data <- create_germline_qval(data, params, sex)

print(get_mutation_stats(data))

# The germline is then filtered using log10(qval) < log10(0.00001), which
# represents a q-value of 0.00001. This is equivalent to a p-value of 0.99999.
# This is the default value used in the original script, but it can be changed
# using the germline_cutoff parameter in the params file.
# This represents a very stringent filter, with only the most confident
# (99.999% certain) mutations being retained. This is to ensure that the
# mutations used to construct the phylogenetic tree are as accurate as
# possible.
print("Filtering germline mutations and mutations of insufficient depth")
# Apply depth and germline filters using params file

data <- data %>%
    mutate(
        removed_by_depth_germline = ifelse(
            (sufficient_depth & log10(germline_qval) < params$germline_cutoff),
            FALSE, TRUE
        )
    )


if (params$beta_binom_shared) {
    print("Flagging shared mutations")
    data <- flag_shared_mutations(data)
}

data <- data %>%
    filter(removed_by_depth_germline == FALSE & shared == TRUE)

print(get_mutation_stats(data))

print("Estimating rho (overdispersion) for each mutation")
data <- estimate_beta_binomial_rho(data, params)

print("Filtering mutations with high rho")

data <- data %>%
    mutate(
        filtered_by_overdispersion = ifelse(
            (Mutation_Type == "SNV" & rho_estimate < params$snv_rho) |
            (Mutation_Type == "indel" & rho_estimate < params$indel_rho),
            FALSE, TRUE
        )
    )


data <- data %>%
    filter(filtered_by_overdispersion == FALSE)

print(get_mutation_stats(data))


# Filter out non-autosomal mutations and those with fewer than three supporting
# reads
data <- data %>%
    mutate(
        filtered_by_autosome_and_support_reads = ifelse(
            (is_XY_chromosomal | NR < 3),
            TRUE, FALSE
        )
    )

data <- data %>%
    filter(filtered_by_autosome_and_support_reads == FALSE)

print(get_mutation_stats(data))

test <- flag_close_to_indel(data, params)

test <- test %>%
    filter(close_to_indel == FALSE)

print(get_mutation_stats(test))

# Testing up to here complete.

# Avoid numerical error from dividing by 0.
data <- data %>%
    # set NR == 0 to 1
    mutate(NR = ifelse(NR == 0, 1, NR))




print("Fitting binomial mixture models")

mix_model_results <- mixture_modelling(data %>% filter(NV > 0), params)

plot_mixture_models(mix_model_results, params)

data <- data %>%
    mutate(
        filtered_by_mixture_model = ifelse(
            (peak_VAF > params$vaf_threshold_mixmodel),
            FALSE, TRUE
        )
    )


stats <- get_mutation_stats(data)
print(stats)


# SAVE FILTERED DATA TO FILE
write.table(
    data,
    file = paste0(
        params$output_dir, "/", params$donor_id, "_final_filtered_data.csv"
    ),
    sep = ",", quote = FALSE, row.names = FALSE
)

# load from disk
data <- read.table(
    file = paste0(
        params$output_dir, "/", params$donor_id, "_final_filtered_data.csv"
    ),
    header = TRUE, sep = ",", stringsAsFactors = FALSE
)

print("Constructing fasta file")

# Add 0 to the NV and the VAF, and 1 to the NR for the mutations not present in one sample
# but present in another. This is to avoid numerical errors when creating the
# genotype.

genotyping_data <- data %>%
    ungroup() %>%
    select(Muts, Sample, NV, NR, VAF) %>%
    complete(Muts, Sample, fill = list(NV = 0, NR = 1, VAF = 0)) %>%
    arrange(Muts, Sample)

stats <- get_mutation_stats(genotyping_data )
print(stats)



almost_binary_genotypes <- create_binary_genotype(
    genotyping_data, sex, params
)
make_fasta(almost_binary_genotypes, params)

print("Running MPBoot")
if (params$only_snvs) {
    command <- paste0(
        params$mpboot_path, " -s ", params$output_dir, "/",
        params$donor_id, "_SNV_for_MPBoot.fa -bb 1000"
    )
} else {
    command <- paste0(
        params$mpboot_path, " -s ", params$output_dir, "/",
        params$donor_id, "_indel_for_MPBoot.fa -bb 1000"
    )
}

system(command, ignore.stdout = T)  # TODO Should we ignore stdout?

mpboot_out <- paste0(
    params$output_dir, "/",
    params$donor_id, "_SNV_for_MPBoot.fa.treefile"
)
tree <- prepare_tree(mpboot_out, params)

assign_mutations_and_plot(
    tree, data, almost_binary_genotypes, params
)

