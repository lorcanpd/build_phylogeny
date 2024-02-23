# read_params.R
library(jsonlite)



# TODO: Add something to deal with when NR and NV matrices are not supplied and
#  cgpvaf_output is used instead as can be a list when default type is
#  character.

# Function to check and enforce types
enforce_types <- function(params, default_params) {

    for (param_name in names(params)) {
        if (param_name %in% names(default_params)) {
            # Check if types match, allowing integers where numerics are expected
            expected_type <- class(default_params[[param_name]])
            actual_type <- class(params[[param_name]])

            # If expected is numeric and actual is integer, or vice versa,
            # consider it okay
            if (!((expected_type == "numeric" && actual_type == "integer") ||
                (expected_type == "integer" && actual_type == "numeric"))) {
                # If not the special case above, and types still don't match,
                # issue warning
                if (actual_type != expected_type) {
                    warning(
                        sprintf(
                            "Type mismatch for parameter '%s'. Expected %s, got %s. Using default value.",
                            param_name, expected_type, actual_type
                        )
                    )
                    # Use the default value in case of mismatch
                    params[[param_name]] <- default_params[[param_name]]
                }
            }
        }
    }
    return(params)
}

# if(is.null(opt$exclude_samples)) {samples_exclude=NULL} else {samples_exclude=unlist(strsplit(x=opt$exclude_samples,split = ","))}
# if(is.null(opt$cnv_samples)) {samples_with_CNVs=NULL} else {samples_with_CNVs=unlist(strsplit(x=opt$cnv_samples,split = ","))}
# if(is.null(opt$cgpvaf_output)) {cgpvaf_paths=NULL} else {cgpvaf_paths=unlist(strsplit(x=opt$cgpvaf_output,split = ","))}


# TODO: Remove some default values? Would we rather pipeline throw an error
#  rather than run?
read_json_params <- function(json_file) {
    # Read parameters from JSON file
    if (file.exists(json_file)) {
        file_params <- jsonlite::fromJSON(json_file)
    } else {
        file_params <- list()
    }

        # Have removed NULL values to that types can be enforced.
    default_params <- list(
        donor_id = "Patient", # Patient/donor ID to add to names of output files
        input_nv = "", # Input NV matrix (rows are variants, columns are samples)
        input_nr = "", # Input NR matrix (rows are variants, columns are samples)
        # CGPVaf output file, instead of NR/NV matrices - can be multiple files,
        # i.e. indel and snv data for the same donor (comma-separated)
        cgpvaf_output = list(),
        output_dir = "", # Output directory for files
        # Only run beta-binomial filter on shared mutations. If FALSE, run on all
        # mutations, before germline/depth filtering
        beta_binom_shared = TRUE,
        ncores = 1, # Number of cores to use for the beta-binomial step
        # Name of the dummy normal to exclude from cgpVAF output
        normal_flt = "PDv37is",
        snv_rho = 0.1, # Rho value threshold for SNVs
        indel_rho = 0.15, # Rho value threshold for indels
        min_cov = 10, # Lower threshold for mean coverage across variant site
        max_cov = 500, # Upper threshold for mean coverage across variant site
        # If indel file is provided, only use SNVs to construct the tree (indels
        # will still be mapped to branches)
        only_snvs = TRUE,
        # If both indels and SNVs are provided, plot trees separately for each
        split_trees = TRUE,
        # Keep an ancestral branch in the phylogeny for mutation mapping
        keep_ancestral = FALSE,
        # Option to manually exclude certain samples from the analysis, separate
        # with a comma
        exclude_samples = list(),
        # Samples with CNVs, exclude from germline/depth-based filtering, separate
        # with a comma
        cnv_samples = list(),
        # VAF threshold (autosomal) below which a variant is absent
        vaf_absent = 0.1,
        # VAF threshold (autosomal) above which a variant is present
        vaf_present = 0.3,
        # Use a binomial mixture model to filter out non-clonal samples?
        mixmodel = FALSE,
        treemut_pval = 0.01, # P-value threshold for treemut's mutation assignment
        # Use a binomial mixture model to filter out non-clonal samples?
        genotype_conv_prob = FALSE,
        # P-value threshold for somatic presence if generating a probabilistic
        # genotype matrix
        min_pval_for_true_somatic = 0.05,
        # Minimum variant reads used in generating a probabilistic genotype matrix
        min_variant_reads_shared = 2,
        # Minimum VAF used in generating a probabilistic genotype matrix
        min_vaf_shared = 2,
        # Convert dichotomous tree from MPBoot to polytomous tree
        create_multi_tree = TRUE,
        mpboot_path = "",  # Path to MPBoot executable
        germline_cutoff = -5, # Log10 of germline qval cutoff
        # Reference genome fasta for plotting mutational spectra
        genome_file = "/nfs/cancer_ref01/Homo_sapiens/37/genome.fa",
        # Plot mutational spectra?
        plot_spectra = FALSE,
        # Maximum number of SNVs to plot in mutational spectra
        max_muts_plot = 5000,
        # VAF threshold for the mixture modelling step to consider a sample clonal
        vaf_treshold_mixmodel = 0.3
    )

    file_params <- enforce_types(file_params, default_params)

    # Combine default parameters with file parameters
    # File parameters override default parameters if they exist
    combined_params <- modifyList(default_params, file_params)

    if (is.null(combined_params$input_nv) & is.null(combined_params$input_nr) & is.null(combined_params$cgpvaf_output)) {
        stop("No input NV/NR matrices or CGPVaf output files provided.")
    }
    # TODO: Make this simpler by overwriting exclude_samples rather than
    #  creating a new variable.
    if (length(combined_params$exclude_samples) == 1 &&
        combined_params$exclude_samples[[1]] == "") {
        combined_params$samples_exclude <- list(character(0))
    } else {
            combined_params$samples_exclude <- combined_params$exclude_samples
    }

    if (length(combined_params$cnv_samples) == 1 &&
        combined_params$cnv_samples[[1]] == "") {
        combined_params$samples_with_CNVs <- list(character(0))
    } else {
            combined_params$samples_with_CNVs <- combined_params$cnv_samples
    }

    if (length(combined_params$cgpvaf_output) == 1 &&
        combined_params$cgpvaf_output[[1]] == "") {
        combined_params$cgpvaf_paths <- list(character(0))
    } else {
            combined_params$cgpvaf_paths <- combined_params$cgpvaf_output
    }

    return(combined_params)
}
