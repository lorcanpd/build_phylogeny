







packages <- c(
    "ggplot2", "ape", "seqinr", "stringr", "data.table", "tidyr",  "dplyr",
    "VGAM", "MASS", "ggtree", "Rsamtools", "GenomicRanges", "tidyverse"
)

for (package in packages) {
    if (!require(package, character.only = TRUE)) {
        library(package, character.only = TRUE)
    }
}



params <- list(
    input_nr = "unit_tests/NR_dummy.txt",
    input_nv = "unit_tests/NV_dummy.txt",
    cgpvaf_paths = "",
    samples_with_CNVs = c(""), # c("Sample_2"),
    min_cov = 10,
    max_cov = 500,
    ncores = 1,
    beta_binom_shared = TRUE,
    germline_cutoff = -5,
    snv_rho = 0.1,
    indel_rho = 0.15,
    vaf_absent = 0.1,
    vaf_present = 0.3,
    vaf_treshold_mixmodel = 0.1,
    # Minimum variant reads used in generating a probabilistic genotype matrix
    min_variant_reads_shared = 2,
    # Minimum VAF used in generating a probabilistic genotype matrix
    min_vaf_shared = 2,
    genotype_conv_prob = TRUE,
    only_snvs = TRUE,
    low_vaf_threshold = 0.1 # SET THIS WITH MANAS - 0.1 for now
)

print("Reading data")
source("read_data.R")
# source("generate_dummy_data.R")
data <- process_data(params)

print("Generating filtering flags")
source("filtering_flags.R")
data <- create_XY_flag(data)
sex <- determine_sex(data)
data <- create_CNV_flag(data, params)
data <- determine_mutation_type(data)
data <- create_depth_flag(data, params, sex)

print("Estimating germline q-values")
source("binomial_testing.R")
data <- create_germline_qval(data, params, sex)

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
    filter(
        # Depth filter
        sufficient_depth == TRUE,
        # Germline filter
        log10(germline_qval) < params$germline_cutoff
    )

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
# Filter out mutations with rho > 0.1 for SNVs and rho > 0.15 for indels
filtered_data <- filtered_data %>%
    filter(
        (Mutation_Type == "SNV" & rho_estimate < params$snv_rho) |
        (Mutation_Type == "indel" & rho_estimate < params$indel_rho)
    )

source("mixture_modelling.R")
# Filter out non-autosomal mutations and those with fewer than three supporting
# reads
mixture_model_data <- filtered_data %>%
    filter(
        !is_XY_chromosomal,
        NR >= 3
    )

print("Fitting binomial mixture models")
mixture_model_results <- fit_binom_mix_model(
    mixture_model_data, params, tolerance = 1e-6, max_iter = 5000
)
print("Finished fitting binomial mixture models")

mixture_model_data <- get_peak_VAFs(mixture_model_data, mixture_model_results)

# plot_all_mixture_models(mixture_model_results, mixture_model_data)

filtered_mixture_model_data <- mixture_model_data %>%
    filter(
        peak_VAF > params$vaf_treshold_mixmodel
    )


print("Constructing fasta file")

source("make_fasta.R")
almost_binary_genotypes <- create_binary_genotype(
    filtered_mixture_model_data, sex, params
)
make_fasta(almost_binary_genotypes, params)


# TODO write plot spectra code which can be applied to filtered data and can be
#  used to plot spectra for the data which has been filtered by any filtering
#  step.
# e.g. filtered_out <- data %>% filter(!sufficient_depth) etc.


print("Running MPBoot")
if (params$only_snvs) {
    command <- paste0(
        params$mpboot_path, "mpboot -s ", params$output_dir, "/",
        params$donor_id, "_SNV_for_MPBoot.fa -bb 1000"
    )

} else {
    command <- paste0(
        params$mpboot_path, "mpboot -s ", params$output_dir, "/",
        params$donor_id, "_indel_for_MPBoot.fa -bb 1000"
    )
}

system(command, ignore.stdout = T)  # TODO Should we ignore stdout?

source("mapping_trees.R")
# TODO WHERE DOES MPBOOT OUTPUT GO?
mpboot_out <- ""
tree <- prepare_tree(mpboot_out, params)

assign_mutations_and_plot(tree, almost_binary_genotypes, params)




