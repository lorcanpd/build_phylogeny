
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
blood_data <- process_data(blood_params)
gut_data <- process_data(gut_params)
data <- bind_rows(blood_data, gut_data)

print("Reading in unfiltered data for read ration comparison")
unfiltered_blood_params <- read_json_params("unfiltered_blood_params.json")
unfiltered_gut_params <- read_json_params("unfiltered_gut_params.json")
unfiltered_blood_data <- process_data(unfiltered_blood_params)
unfiltered_gut_data <- process_data(unfiltered_gut_params)
unfiltered_data <- bind_rows(unfiltered_blood_data, unfiltered_gut_data)

# Reduce memory usage by removing unneeded data
rm(
    blood_data,
    gut_data,
    unfiltered_blood_data,
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


data <- flag_binary_germline(data, presence_threshold = 3)

test <- flag_binary_shared(data, presence_threshold = 3)
test <- flag_binary_unique(test, presence_threshold = 3)

test %>%
    group_by(Sample) %>%
    summarise(
        total = n(),
        below_thresholds = sum(!binary_unique & !binary_shared),
        unique = sum(binary_unique),
        shared = sum(binary_shared)
    ) %>%
    select(Sample, total, below_thresholds, shared, unique) %>%
    print()

data <- data %>%
    group_by(Muts) %>%
    filter(all(binary_germline) == FALSE & all(sufficient_depth)) %>%
    ungroup()

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

data <- flag_close_to_indel(data, blood_params)

data <- data %>%
    filter(close_to_indel == FALSE)

data <- flag_binary_shared(data, presence_threshold = 3)
data <- flag_binary_unique(data, presence_threshold = 3)

data %>%
    group_by(Sample) %>%
    summarise(
        total = n(),
        below_thresholds = sum(!binary_unique & !binary_shared),
        unique = sum(binary_unique),
        shared = sum(binary_shared)
    ) %>%
    select(Sample, total, below_thresholds, shared, unique) %>%
    print()

# Avoid numerical error from dividing by 0.
data <- data %>%
    mutate(NR = ifelse(NR == 0, 1, NR))


print("Fitting binomial mixture models")
# Test how changing the addition of 0 NVs effects the behaviour of the mixture
# model. In particular how does this alter the output?
data <- mixture_modelling(data %>% filter(NV > 0), blood_params)

plot_mixture_models(data, blood_params)

data <- data %>%
    mutate(
        filtered_by_mixture_model = ifelse(
            (peak_VAF > blood_params$vaf_threshold_mixmodel),
            FALSE, TRUE
        )
    )

# IMPORTANT!!!
# Ignore this for now as this would remove all the samples haha. Once all
# samples are included then uncomment.
# data <- data %>%
#     filter(filtered_by_mixture_model == FALSE)

print(get_mutation_stats(data))


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
    select(Muts, Sample, NV, NR, VAF, Mutation_Type, is_XY_chromosomal) %>%
    complete(
        Muts, Sample,
        fill = list(
            NV = 0, NR = 1, VAF = 0, Mutation_Type = "SNV",
            is_XY_chromosomal = FALSE
        )
    ) %>%
    arrange(Muts, Sample)

print(get_mutation_stats(genotyping_data))

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

system(command, ignore.stdout = F)

mpboot_out <- paste0(
    blood_params$output_dir, "/",
    blood_params$donor_id, "_SNV_for_MPBoot.fa.treefile"
)
tree <- prepare_tree(mpboot_out, blood_params)

# Debug this.
assign_mutations_and_plot(
    tree, genotyping_data, almost_binary_genotypes, blood_params
)

