
packages <- c(
    "data.table", "dplyr", "tidyr", "magrittr", "readr", "stringr", "tibble"
)

for (package in packages) {
    if (!require(package, character.only = TRUE)) {
        library(package, character.only = TRUE)
    }
}

# process files input data, taking in a cgpvaf file or a set of NV/NR matrices
# and outputting a list of dataframes with the following columns:
# Sample, Chr, Pos, Ref, Alt, NV, NR, VAF, Gene, Depth, Impact, AAchange, Snp, Flag, ID

process_cgpvaf_file <- function(params) {
    # Function to process a single cgpvaf output file. This function is used
    # when only a single file is provided in the params file. The function
    # returns a list containing the Muts, NR, and NV objects.
    data <- fread(params$cgpvaf_paths, header = TRUE, data.table = FALSE)
    Muts <- paste(data$Chrom, data$Pos, data$Ref, data$Alt, sep = "_")
    NR <- data[, grepl("DEP", colnames(data)) & !grepl(paste(
        c(params$normal_flt, params$samples_exclude), collapse = "|"
    ), colnames(data)), with = FALSE]
    NV <- data[, grepl("MTR", colnames(data)) & !grepl(paste(
        c(params$normal_flt, params$samples_exclude), collapse = "|"
    ), colnames(data)), with = FALSE]

    list(Muts = Muts, NR = NR, NV = NV)
}

process_multiple_cgpvaf_files <- function(params) {
    results <- lapply(
        params$cgpvaf_paths,
        process_cgpvaf_file,
        params$normal_flt,
        params$samples_exclude
    )
    Muts_combined <- unlist(lapply(results, '[[', 'Muts'))
    NR_combined <- do.call(rbind, lapply(results, '[[', 'NR'))
    NV_combined <- do.call(rbind, lapply(results, '[[', 'NV'))
    list(Muts = Muts_combined, NR = NR_combined, NV = NV_combined)
}

process_nr_nv_files <- function(params) {
    NR <- fread(params$input_nr, data.table = F)
    rownames(NR) <- NR[, 1]
    NR <- NR[, -1]

    NV <- fread(params$input_nv, data.table = F)
    rownames(NV) <- NV[, 1]
    NV <- NV[, -1]

    list(Muts = rownames(NV), NR = NR, NV = NV)
}


convert_to_long_table <- function (data) {
    # Reshape the data into a long format table
    data_long <- data$NV %>%
        rownames_to_column("Muts") %>%
        gather(Sample, NV, -Muts) %>%
        left_join(
            data$NR %>%
                rownames_to_column("Muts") %>%
                gather(Sample, NR, -Muts),
            by = c("Muts", "Sample")
        ) %>%
        left_join(
            data$Muts %>%
                as.data.frame() %>%
                setNames("Muts"),  # TODO RENAME TO MUT_ID?
            by = "Muts"
        ) %>%
        mutate(
            Chr = sapply(strsplit(Muts, "_"), `[`, 1),
            Pos = as.numeric(sapply(strsplit(Muts, "_"), `[`, 2)),
            Ref = sapply(strsplit(Muts, "_"), `[`, 3),
            Alt = sapply(strsplit(Muts, "_"), `[`, 4),
            VAF = ifelse(NR == 0, 0, NV / NR)
        )
    return(data_long)
}


process_data <- function(params) {
    if (nzchar(params$cgpvaf_paths)) {
        if (length(params$cgpvaf_paths) == 1) {
            result <- process_cgpvaf_file(params)
        } else {
            result <- process_multiple_cgpvaf_files(params)
        }
    } else {
        if (nzchar(params$input_nr) & nzchar(params$input_nv)) {
            result <- process_nr_nv_files(params)
        } else {
            stop("Please provide either NV and NR files or a path to CGPVaf output")
        }
    }

    # Set the rownames and colnames for NR and NV
    rownames(result$NV) <- rownames(result$NR) <- result$Muts
    samples <- colnames(result$NR) <- colnames(result$NV) <- gsub(
        "_DEP", "", colnames(result$NR)
    )

    data <- list(
        Muts = result$Muts, NR = result$NR, NV = result$NV, Sample = samples
    )
    data_long <- convert_to_long_table(data)

    return(data_long)
}
