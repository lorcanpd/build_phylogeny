# library(tidyverse)
library(ggplot2)
library(gridExtra)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
# Call in package which allows scanFa
library(Rsamtools)

annotate_mutations <- function(mutations, genomeFile) {

    ntcomp <- c(T = "A", G = "C", C = "G", A = "T")

    mutations_accross_samples <- mutations %>%
        select(Muts, Chr, Pos, Ref, Alt) %>%
        distinct()

    process_mutation_chunk <-function(mutation_chunk, genomeFile) {
        mutation_chunk %>%
            mutate(
                trinuc_ref = as.vector(
                    scanFa(
                        genomeFile, GRanges(Chr, IRanges(Pos - 1, Pos + 1))
                    )
                ),
                sub = if_else(
                    Ref %in% c("A", "G"),
                    paste(ntcomp[Ref], ntcomp[Alt], sep = ">"),
                    paste(Ref, Alt, sep = ">")
                ),
                trinuc_ref_py = if_else(
                    Ref %in% c("A", "G"),
                    stringr::str_to_upper(
                        sapply(
                            strsplit(trinuc_ref, ""),
                            function(x) paste(rev(x), collapse = "")
                        )
                    ) %>%
                        stringr::str_replace_all(ntcomp),
                    trinuc_ref
                )
            )
    }

    num_cores <- ifelse(is.null(params$ncores), 1, params$ncores)

    if (num_cores > 1) {
        library(parallel)
        library(doParallel)
        cl <- makeCluster(params$ncores)
        registerDoParallel(cl)
        # Split data into chunks
        mutations_chunks <- split(
            mutations_accross_samples,
            cut(seq(nrow(mutations_accross_samples)), num_cores)
        )
        # Use mclapply for parallel processing
        results <- mclapply(
            mutations_chunks, process_mutation_chunk, genomeFile,
            mc.cores = num_cores
        )
        # Combine results
        annotated_mutations <- do.call(rbind, results)
        stopCluster(cl)
    } else {
        # If only one core, just apply the function normally
        annotated_mutations <- process_mutation_chunk(
            mutations_accross_samples, genomeFile
        )
    }

    # join annotated_mutations to original data
    mutations <- annotated_mutations %>%
        left_join(
            mutations,
            by = c("Muts", "Chr", "Pos", "Ref", "Alt")
        )

    return(mutations)
}

get_ordered_trinucleotides <- function() {
    # Mutation types
    sub_types <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

    # Bases
    bases <- c("A", "C", "G", "T")

    # Generate all combinations of flanking bases
    contexts <- expand.grid(
        base1 = bases, base2 = bases, stringsAsFactors = FALSE
    )

    trinucleotides <- c()

    # Loop through contexts and for each context, append all subtypes
    for (i in 1:nrow(contexts)) {
        for (sub in sub_types) {
            ctx <- with(contexts, paste0(base1[i], substr(sub, 1, 1), base2[i]))
            trinucleotides <- c(trinucleotides, paste0(ctx, ",", sub))
        }
    }

    return(trinucleotides)
}

count_substitutions <- function(mutations) {

    ordered_trinucleotides <- get_ordered_trinucleotides()

    freqs <- mutations %>%
        group_by(trinuc_ref_py, sub) %>%
        summarise(freq = n(), .groups = "drop") %>%
        mutate(full=paste0(trinuc_ref_py, ",", sub)) %>%
        complete(full = ordered_trinucleotides, fill = list(freq = 0)) %>%
        filter(full %in% ordered_trinucleotides) %>%
        select(
            full, trinucleotide = trinuc_ref_py, substitution = sub, freq
        ) %>%
        mutate(
            trinucleotide = if_else(
                is.na(trinucleotide), substr(full, 1, 3), trinucleotide
            ),
            substitution = if_else(
                is.na(substitution), substr(full, 5, 7), substitution
            )
        ) %>%
        arrange(substitution, trinucleotide) %>%
        mutate(full = factor(full, levels = unique(full)))

    return(freqs)
}


plot_data <- function(freqs, sample_name, add_to_title = "") {

    colvec <- c("C>A" = "dodgerblue", "C>G" = "black", "C>T" = "red",
            "T>A" = "grey70", "T>C" = "olivedrab3", "T>G" = "plum2")
    sub_vec <- names(colvec)
    # Round max(freqs$freq) to the nearest 100
    top <- round(max(freqs$freq) * 1.25, -2)

    # Create the plot
    plot <- ggplot(freqs, aes(x = full, y = freq, fill = substitution)) +
        geom_bar(stat = "identity", show.legend = FALSE, width = 0.65) +
        scale_fill_manual(values = colvec) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.spacing.x = unit(0.5, "lines"),
            axis.line.y = element_line(color = "black"),
            axis.ticks.y = element_line(color = "black")
        ) +
        labs(
            y = "Number of mutations",
            x = "",
            title = paste0("Number of mutations: ", sum(freqs$freq))
        ) +
        scale_x_discrete(
            labels = function(x) {
                ifelse(
                    seq_along(x) %% 2 == 0, "", substr(x, 1, 3)
                )
            }
        ) +
        scale_y_continuous(
            expand = c(0, 0),
            breaks = seq(0, top, top%/%5)
        ) +
        coord_cartesian(xlim = c(-1, NA))

        # Calculate label positions
    num_trinuc_per_sub <- length(unique(freqs$full)) / length(sub_vec)
    label_positions <- seq(
        (num_trinuc_per_sub + 0.5) / 2,
        by = num_trinuc_per_sub,
        length.out = length(sub_vec)
    )

    # Add rectangles and labels using annotate
    for (i in seq_along(sub_vec)) {
        xpos1 <- label_positions[
            i] - num_trinuc_per_sub / 2 + num_trinuc_per_sub * 0.025
        xpos2 <- label_positions[
            i] + num_trinuc_per_sub / 2 - num_trinuc_per_sub * 0.025
        ypos1 <- max(freqs$freq) * 1.05
        ypos2 <- ypos1 + max(freqs$freq) * 0.1

        plot <- plot +
            annotate(
                "rect", xmin = xpos1, xmax = xpos2, ymin = ypos1, ymax = ypos2,
                fill = colvec[sub_vec[i]], color = NA
            ) +
            annotate(
                "text", x = mean(c(xpos1, xpos2)),
                # y = ypos1 + max(freqs$freq) * 0.05,
                y = ypos1 + 1.5 * (ypos2 - ypos1),
                label = sub_vec[i], size = 3
            ) +
            expand_limits(y = top)
    }

    return(plot)
}

plot_spectrum_lattice <- function(
    data, save, add_to_title = "", by_sample = FALSE
) {
    if (by_sample) {
        plots <- data %>%
            group_by(Sample) %>%
            do({
                freqs_full <- count_substitutions(.)
                plot_data(freqs_full, unique(sample_data$Sample), add_to_title)
            })
        plot_name <- paste0("plots/spectrum_by_sample", add_to_title, ".pdf")
    } else {
        plots <- list()
        freqs_full <- count_substitutions(data)

        plots[[1]] <- plot_data(
            freqs_full, "All samples", add_to_title
        )
        plot_name <- paste0("plots/all_sample_spectra", add_to_title, ".pdf")
    }

    lattice_plot <- wrap_plots(plots, ncol = 1)  # Adjust ncol as needed

    if (save) {
        ggsave(plot_name, lattice_plot, width = 12, height = 4 * length(plots))
    } else {
        print(lattice_plot)
    }
}



plot_vaf_over_chromosomes <- function (data, sex, sample_name) {

    chromosomes <- paste0("chr", 1:22)
    if (sex == "male") {
        # Adjust the frequency of chrY based on the y_ratio
        chromosomes <- c(chromosomes, "chrX", "chrY")
    } else {
        chromosomes <- c(chromosomes, "chrX")
    }

    data <- data %>%
        filter(Sample == sample_name) %>%
        mutate(Chr = factor(Chr, levels = chromosomes)) %>%
        arrange(Chr, Pos)

    plots <- list()

    # Calculate total range (or total number of data points) for normalization
    total_range <- nrow(data)
    widths <- c()
    # Create a plot for each chromosome and store in the list
    for (chr in levels(data$Chr)) {
        chr_data <- data[data$Chr == chr,]
        p <- ggplot(chr_data, aes(x = Pos, y = VAF)) +
            ylim(0, 1) +
            geom_point(size=0.2) +
            ggtitle(chr)

        # Apply different theme settings for the first chromosome
        if (chr == levels(data$Chr)[1]) {
            p <- p + theme_minimal()
        } else {
            p <- p +
                theme_minimal() +
                theme(
                    axis.title.y = element_blank(),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank()
                )
        }

        plots[[chr]] <- p

        # Calculate width based on range
        chr_range <- nrow(chr_data) #max(chr_data$Pos) - min(chr_data$Pos)
        widths <- c(widths, chr_range / total_range)
    }

    # Combine plots with calculated widths
    do.call(grid.arrange, c(plots, ncol = 1, list(widths = widths)))

    return(plots)
}



