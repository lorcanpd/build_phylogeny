
# This script showcases how the modularised pipeline can be used to build a
# phylogenetic tree, with optional iterative filtering and plotting of
# mutational signatures. The pipeline is designed to be flexible and modular,
# allowing for easy modification and extension. Uses are expected to write a
# custom pipeline which uses the functionality demonstrated here to suit their
# specific needs.

library(buildphylogeny)

print("Reading in mapping and base quality filtered data")
blood_params <- read_json_params("blood_params.json")
gut_params <- read_json_params("gut_params.json")
# blood_data <- process_data(blood_params)
gut_data <- process_data(gut_params)
# data <- bind_rows(blood_data, gut_data)
data <- gut_data

print("Reading in unfiltered data for read ration comparison")
# unfiltered_blood_params <- read_json_params("unfiltered_blood_params.json")
unfiltered_gut_params <- read_json_params("unfiltered_gut_params.json")
# unfiltered_blood_data <- process_data(unfiltered_blood_params)
unfiltered_gut_data <- process_data(unfiltered_gut_params)
# unfiltered_data <- bind_rows(unfiltered_blood_data, unfiltered_gut_data)
unfiltered_data <- unfiltered_gut_data

rm(
    # blood_data,
    gut_data,
    # unfiltered_blood_data,
    unfiltered_gut_data
)

print("Using read ratio flags to filter out mutations")
data <- read_ratio_flags(data, unfiltered_data, ratio_threshold = 0.75)

data <- data %>%
    filter(filtered_by_read_ratio == FALSE)


data <- annotate_mutations(data, blood_params)

print("Generating filtering flags")
data <- create_XY_flag(data)
sex <- determine_sex(data)
data <- create_CNV_flag(data, blood_params)
data <- determine_mutation_type(data)
data <- create_depth_flag(data, blood_params, sex)

print("Estimating germline q-values")
data <- create_germline_qval(data, blood_params, sex)

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
            (sufficient_depth &
                log10(germline_qval) < blood_params$germline_cutoff),
            FALSE, TRUE
        )
    )


if (blood_params$beta_binom_shared) {
    print("Flagging shared mutations")
    data <- flag_shared_mutations(data)
}

data <- data %>%
    filter(removed_by_depth_germline == FALSE & shared == TRUE)

print(get_mutation_stats(data))

print("Estimating rho (overdispersion) for each mutation")
data <- estimate_beta_binomial_rho(data, blood_params)

print("Filtering mutations with high rho")

data <- data %>%
    mutate(
        filtered_by_overdispersion = ifelse(
            (Mutation_Type == "SNV" & rho_estimate < blood_params$snv_rho) |
            (Mutation_Type == "indel" & rho_estimate < blood_params$indel_rho),
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

data <- flag_close_to_indel(data, blood_params)

data <- data %>%
    filter(close_to_indel == FALSE)

print(get_mutation_stats(test))

# Avoid numerical error from dividing by 0.
data <- data %>%
    mutate(NR = ifelse(NR == 0, 1, NR))




print("Fitting binomial mixture models")

# Test how changing the addition of 0 NVs effects the behaviour of the mixture
# model. In particular how does this alter the output?
mix_model_results <- mixture_modelling(data %>% filter(NV > 0), blood_params)

plot_mixture_models(mix_model_results, blood_params)

data <- data %>%
    mutate(
        filtered_by_mixture_model = ifelse(
            (peak_VAF > blood_params$vaf_threshold_mixmodel),
            FALSE, TRUE
        )
    )


stats <- get_mutation_stats(data)
print(stats)


# SAVE FILTERED DATA TO FILE
write.table(
    data,
    file = paste0(
        blood_params$output_dir, "/",
        blood_params$donor_id, "_final_filtered_data.csv"
    ),
    sep = ",", quote = FALSE, row.names = FALSE
)

# load from disk
data <- read.table(
    file = paste0(
        blood_params$output_dir, "/",
        blood_params$donor_id, "_final_filtered_data.csv"
    ),
    header = TRUE, sep = ",", stringsAsFactors = FALSE
)

print("Constructing fasta file")

if (blood_params$genotype_conv_prob) {
    data <- flag_conv_shared_mutations(data, blood_params)
    # TODO: Find out what Tim meant in his original script with min_shared_vaf
    #  default value = 2 and then indexing this integer (which will fail) to
    #  get separate values for autosomal and XY depending on sample sex. In
    #  any case we are not using it here.
}

genotyping_data <- data %>%
    ungroup() %>%
    select(Muts, Sample, NV, NR, VAF, Mutation_Type) %>%
    complete(
        Muts, Sample,
        fill = list(NV = 0, NR = 1, VAF = 0, Mutation_Type = "SNV")
    ) %>%
    arrange(Muts, Sample)

stats <- get_mutation_stats(genotyping_data)
print(stats)

almost_binary_genotypes <- create_binary_genotype(
    genotyping_data, sex, blood_params
)
make_fasta(almost_binary_genotypes, blood_params)

print("Running MPBoot")
if (blood_params$only_snvs) {
    command <- paste0(
        blood_params$mpboot_path, " -s ", blood_params$output_dir, "/",
        blood_params$donor_id, "_SNV_for_MPBoot.fa -bb 1000"
    )
} else {
    command <- paste0(
        blood_params$mpboot_path, " -s ", blood_params$output_dir, "/",
        blood_params$donor_id, "_indel_for_MPBoot.fa -bb 1000"
    )
}

system(command, ignore.stdout = T)  # TODO Should we ignore stdout?

mpboot_out <- paste0(
    blood_params$output_dir, "/",
    blood_params$donor_id, "_SNV_for_MPBoot.fa.treefile"
)
tree <- prepare_tree(mpboot_out, blood_params)

# Debug this.
assign_mutations_and_plot(
    tree, data, almost_binary_genotypes, blood_params
)

