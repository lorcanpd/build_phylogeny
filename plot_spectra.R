library(tidyverse)
library(GenomicRanges)
library(ggplot2)
library(patchwork)

annotate_mutations <- function(mutations, genomeFile) {
    ntcomp <- c(T = "A", G = "C", C = "G", A = "T")
    mutations %>%
        mutate(
            trinuc_ref = scanFa(
                genomeFile, GRanges(Chr, IRanges(Pos - 1, Pos + 1))
            ),
            sub = paste(ref, mut, sep = ">"),
            trinuc_ref_py = if_else(
                ref %in% c("A", "G"),
                stringr::str_replace_all(
                    stringr::str_to_upper(
                        stringr::str_reverse(
                            trinuc_ref
                        )
                    ),
                    ntcomp
                ),
                trinuc_ref
            )
        )
    return(mutations)
    }

count_substitutions <- function(mutations) {
    freqs <- table(
        paste(
            mutations$sub,
            paste(
                substr(mutations$trinuc_ref_py, 1, 1),
                substr(mutations$trinuc_ref_py, 3, 3),
                sep = "-"
            ),
            sep = ","
        )
    )
    sub_vec <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    ctx_vec <- paste(
        rep(c("A", "C", "G", "T"), each = 4),
        rep(c("A", "C", "G", "T"), times = 4),
        sep = "-"
    )
    full_vec <- paste(
        rep(sub_vec, each = 16), rep(ctx_vec, times = 6), sep = ","
    )
    freqs_full <- freqs[full_vec]
    freqs_full[is.na(freqs_full)] <- 0
    names(freqs_full) <- full_vec
    return(freqs_full)
}

plot_data <- function(freqs_full, sample_name, add_to_title = "") {
    xstr <- paste(
        substr(names(freqs_full), 5, 5),
        substr(names(freqs_full), 1, 1),
        substr(names(freqs_full), 7, 7),
        sep = ""
    )
    colvec <- rep(
        c("dodgerblue", "black", "red", "grey70", "olivedrab3", "plum2"),
        each = 16
    )

    plot <- ggplot(
            data.frame(xstr, freqs_full),
            aes(x = xstr, y = freqs_full, fill = xstr)
        ) +
        geom_bar(stat = "identity", color = "black") +
        scale_fill_manual(values = colvec) +
        theme_minimal() +
        labs(
            title = paste0("Sample: ", sample_name, add_to_title),
            x = "",
            y = "Number of mutations"
        ) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    return(plot)
}

plot_spectrum_lattice <- function(
    data, save, genomeFile, add_to_title = "", by_sample = FALSE
) {
    if (by_sample) {
        plots <- data %>%
            group_by(Sample) %>%
            do({
                sample_data <- prepare_data(.)
                annotated_mutations <- annotate_mutations(
                    sample_data, genomeFile
                )
                freqs_full <- count_substitutions(annotated_mutations)
                plot_data(freqs_full, unique(sample_data$Sample), add_to_title)
            })
        plot_name = paste0("spectrum_by_sample", add_to_title, ".pdf")
    } else {
        plots <- list()
        sample_data <- prepare_data(data)
        annotated_mutations <- annotate_mutations(sample_data, genomeFile)
        freqs_full <- count_substitutions(annotated_mutations)
        plots[[1]] <- plot_data(
            freqs_full, "All samples", add_to_title
        )
        plot_name = paste0("all_sample_spectra", add_to_title, ".pdf")
    }

    lattice_plot <- wrap_plots(plots, ncol = 1)  # Adjust ncol as needed

    if (save) {
        ggsave(plot_name, lattice_plot, width = 12, height = 4 * nrow(plots))
    } else {
        print(lattice_plot)
    }
}

