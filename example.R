

# This script showcases how the modularised pipeline can be used to build a
# phylogenetic tree, with optional iterative filtering and plotting of
# mutational signatures. The pipeline is designed to be flexible and modular,
# allowing for easy modification and extension. Uses are expected to write a
# custom pipeline which uses the functionality demonstrated here to suit their
# specific needs.

library(buildphylogeny)
library(treemut)

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

data %>%
    group_by(Sample) %>%
    summarise(
        num_mutations = n()
    ) %>%
    print()

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
filtered_data <- data %>%
    # Group by mutation identifier
    group_by(Muts) %>%
    # Filter groups where at least one sample in the group passes the filters
    filter(
        any(
            sufficient_depth == TRUE &
                log10(germline_qval) < params$germline_cutoff
        )
    ) %>%
    # Ungroup to remove the grouping structure
    ungroup()

filtered_out <- data %>%
    # Group by mutation identifier
    group_by(Muts) %>%
    # Filter to retain mutations if all samples fail the conditions
    filter(
        all(
            sufficient_depth == FALSE |
                log10(germline_qval) >= params$germline_cutoff
        )
    ) %>%
    # Ungroup to remove the grouping structure
    ungroup()

filtered_data %>%
    group_by(Sample) %>%
    summarise(
        num_mutations = n()
    ) %>%
    print()

plot_spectrum_lattice(
    filtered_data, TRUE,
    add_to_title = "_retained_by_depth_germline_filters",
    by_sample = FALSE
)
plot_spectrum_lattice(
    filtered_out, TRUE,
    add_to_title = "_excluded_by_depth_germline_filters",
    by_sample = FALSE
)
rm(filtered_out)

if (params$beta_binom_shared) {
    print("Flagging shared mutations")
    filtered_data <- flag_shared_mutations(filtered_data)
    # Estimate rho for each mutation
    filtered_data <- filtered_data %>%
        filter(shared == TRUE)
}

print("Estimating rho (overdispersion) for each mutation")
filtered_data <- estimate_beta_binomial_rho(filtered_data, params)

print("Filtering mutations with high rho")

filtered_out <- filtered_data %>%
    group_by(Muts) %>%
    filter(
        all(
            (Mutation_Type == "SNV" & rho_estimate >= params$snv_rho) |
            (Mutation_Type == "indel" & rho_estimate >= params$indel_rho)
        )
    ) %>%
    ungroup()

# Filter out mutations with rho > 0.1 for SNVs and rho > 0.15 for indels
filtered_data <- filtered_data %>%
    group_by(Muts) %>%
    filter(
        any(
            (Mutation_Type == "SNV" & rho_estimate < params$snv_rho) |
            (Mutation_Type == "indel" & rho_estimate < params$indel_rho)
        )
    ) %>%
    ungroup()

filtered_data %>%
    group_by(Sample) %>%
    summarise(
        num_mutations = n()
    ) %>%
    print()

plot_spectrum_lattice(
    filtered_data, TRUE,
    add_to_title = "_retained_by_overdispersion_filter",
    by_sample = FALSE
)

plot_spectrum_lattice(
    filtered_out, TRUE,
    add_to_title = "_excluded_by_overdispersion_filter",
    by_sample = FALSE
)
rm(filtered_out)

# Filter out non-autosomal mutations and those with fewer than three supporting
# reads
filtered_data <- filtered_data %>%
    group_by(Muts) %>%
    filter(any(!is_XY_chromosomal & NR >= 3)) %>%
    ungroup()

filtered_out <- filtered_data %>%
    group_by(Muts) %>%
    filter(all(is_XY_chromosomal | NR < 3)) %>%
    ungroup()

filtered_data %>%
    group_by(Sample) %>%
    summarise(
        num_mutations = n()
    ) %>%
    print()

plot_spectrum_lattice(
    filtered_data, TRUE,
    add_to_title = "_retained_by_autosomal_filter",
    by_sample = FALSE
)

plot_spectrum_lattice(
    filtered_out, TRUE,
    add_to_title = "_excluded_by_autosomal_filter",
    by_sample = FALSE
)
rm(filtered_out)

# Avoid numerical error from dividing by 0.
# mixture_model_data <- filtered_data %>%
#     # set NR == 0 to 1
#     mutate(NR = ifelse(NR == 0, 1, NR))


print("Fitting binomial mixture models")
# Currently mixture modellnig can't be done with multiple cores.
# Probably due to the fact I can only get a notebook with 16 gb RAM.
params_for_mixmodel <- params
params_for_mixmodel$ncores <- 1
# mixture_model_data, params, tolerance = 1e-6, max_iter = 5000 # Use these parameters.
mixture_model_results <- fit_binom_mix_model(
    mixture_model_data, params_for_mixmodel, tolerance = 1e-4, max_iter = 500 # Test perameters
)
print("Finished fitting binomial mixture models")

mixture_model_data <- get_peak_VAFs(mixture_model_data, mixture_model_results)

# plot_all_mixture_models(mixture_model_results, mixture_model_data)

filtered_mixture_model_data <- mixture_model_data %>%
    group_by(Muts) %>%
    filter(any(peak_VAF > params$vaf_threshold_mixmodel)) %>%
    ungroup()

filtered_out <- mixture_model_data %>%
    group_by(Muts) %>%
    filter(any(peak_VAF <= params$vaf_threshold_mixmodel)) %>%
    ungroup()


plot_spectrum_lattice(
    filtered_mixture_model_data, TRUE,
    add_to_title = "_retained_by_mixture_model_filter",
    by_sample = FALSE
)

plot_spectrum_lattice(
    filtered_out, TRUE,
    add_to_title = "_excluded_by_mixture_model_filter",
    by_sample = FALSE
)
rm(filtered_out)

# count the number of mutations per sample
filtered_mixture_model_data %>%
    group_by(Sample) %>%
    summarise(
        num_mutations = n()
    ) %>%
    print()

# SAVE FILTERED DATA TO FILE
write.table(
    filtered_mixture_model_data,
    file = paste0(
        params$output_dir, "/", params$donor_id, "_final_filtered_data.csv"
    ),
    sep = ",", quote = FALSE, row.names = FALSE
)

# load from disk
filtered_mixture_model_data <- read.table(
    file = paste0(
        params$output_dir, "/", params$donor_id, "_final_filtered_data.csv"
    ),
    header = TRUE, sep = ",", stringsAsFactors = FALSE
)

print("Constructing fasta file")
almost_binary_genotypes <- create_binary_genotype(
    filtered_mixture_model_data, sex, params
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
    tree, filtered_mixture_model_data, almost_binary_genotypes, params
)

