library(treemut)
library(tidyr)
library(dplyr)
library(tibble)
library(ggtree)


add_ancestral_outgroup <- function(tree, outgroup_name = "Ancestral") {
    # This function adds the ancestral tip at the end
    tmp <- tree$edge
    N <- length(tree$tip.label)
    newroot <- N + 2
    renamedroot <- N + 3
    ancestral_tip <- N + 1
    tmp <- ifelse(tmp > N, tmp + 2, tmp)

    tree$edge <- rbind(c(newroot, renamedroot), tmp, c(newroot, ancestral_tip))
    tree$edge.length <- c(0, tree$edge.length, 0)

    tree$tip.label <- c(tree$tip.label, outgroup_name)
    tree$Nnode <- tree$Nnode + 1
    mode(tree$Nnode) <- "integer"
    mode(tree$edge) <- "integer"
    return(tree)
}


prepare_tree <- function(tree_path, params) {
    tree <- read.tree(tree_path)
    tree <- drop.tip(tree, "Ancestral")

    if (params$keep_ancestral) {
        tree <- add_ancestral_outgroup(tree)
    }

    tree$edge.length <- rep(1, nrow(tree$edge))

    return(tree)
}


adjust_edge_lengths <- function(tree, pval_threshold) {
    # Extract the edges with significant mutation events
    likely_nowhere_else <- table(
        tree$summary$edge_ml[tree$summary$p_else_where < pval_threshold]
    )
    # Check if the tree is an adjusted tree or not and then assign edge lengths
    # and names accordingly.
    if (!is.null(tree$tree$edge)) {
        edge_length <- rep(0, nrow(tree$tree$edge))
        names(edge_length) <- seq_len(nrow(tree$tree$edge))
    } else {
        edge_length <- rep(0, nrow(tree$edge))
        names(edge_length) <- seq_len(nrow(tree$edge))
    }

    # Update edge lengths for edges with significant mutation events
    edge_length[names(likely_nowhere_else)] <- likely_nowhere_else

    # Update the tree object with new edge lengths
    if (!is.null(tree$tree$edge)){
        tree$tree$edge.length <- as.numeric(edge_length)
    } else {
        tree$edge.length <- as.numeric(edge_length)
    }

    return(tree)
}


get_mutations <- function(data, binary_genotypes) {
    # Get mutations from the data
        present_mutations <- binary_genotypes %>%
        group_by(Muts) %>%
        summarise(
            globally_present = sum(binary_genotype > 0) > 0
        ) %>%
        left_join(data, by = "Muts") %>%
        filter(globally_present == TRUE) %>%
        select(Muts, Sample, NV, NR)
    return(present_mutations)
}


reassign_tree <- function(tree, present_mutations, params) {
    tree <- drop.tip(tree, "Ancestral")

    tree$edge.length=rep(1, nrow(tree$edge))
    # genotype_summary <- reconstruct_genotype_summary(tree)

    NR <- present_mutations %>%
        select(Muts, Sample, NR) %>%
        spread(Sample, NR) %>%
        column_to_rownames("Muts") %>%
        as.matrix()
    NV <- present_mutations %>%
        select(Muts, Sample, NV) %>%
        spread(Sample, NV) %>%
        column_to_rownames("Muts") %>%
        as.matrix()

    if (params$keep_ancestral) {
        NR <- cbind(NR, Ancestral = rep(30, nrow(NR)))
        NV <- cbind(NV, Ancestral = rep(0, nrow(NV)))
        error_probabilies <- c(1e-6, rep(0.01, ncol(NV)-1))
        reassigned_tree <- assign_to_tree(
            tree=tree, mtr=NV, dep=NR, error_rate = error_probabilies
        )
    } else {
        reassigned_tree <- assign_to_tree(
            tree=tree,
            mtr=NV, dep=NR
        )
    }
    write_mutations_per_branch(reassigned_tree, NR, params)
    return(reassigned_tree)
}


write_and_plot_tree <- function(tree, params, plot_name) {
    tree_plot <- ggtree(tree$tree) +
        geom_tiplab(aes(x = branch), vjust=-0.3) +
        theme_tree2() +
        xlim(0, max(tree$df$df$expected_edge_length) * 1.3)
    pdf(
        paste0(
            params$output_dir, "/",
            params$donor_id, "_",
            plot_name, "_tree_with_branch_length.pdf"
        )
    )
    print(tree_plot)
    dev.off()

    write.tree(tree$tree, paste0(
        params$output_dir, "/",
        params$donor_id, "_",
        plot_name, "_tree_with_branch_length.tree"
    ))
}


write_mutations_per_branch <- function(tree, NR, params) {
    Mutations_per_branch <- as.data.frame(
        matrix(
            ncol = 4, unlist(strsplit(rownames(NR), split = "_")), byrow = TRUE
        )
    )
    colnames(Mutations_per_branch) <- c("Chr", "Pos", "Ref", "Alt")
    Mutations_per_branch$Branch <- tree$tree$edge[tree$summary$edge_ml, 2]
    Mutations_per_branch <- Mutations_per_branch[
        tree$summary$p_else_where < params$treemut_pval,
    ]

    Mutations_per_branch$Patient <- params$donor_id
    Mutations_per_branch$SampleID <- paste(
        params$donor_id, Mutations_per_branch$Branch, sep = "_"
    )

    if (params$only_snvs) {
        mut_id <- "SNV"
    } else {
        mut_id <- "indel"
    }
    write.table(
        Mutations_per_branch,
        paste0(
            params$output_dir, params$donor_id, "_", mut_id,
            "_assigned_to_branches.tsv"
        ),
        quote = FALSE, row.names = FALSE, sep = "\t"
    )
}


assign_mutations_and_plot <- function (
    tree, data, binary_genotypes, params
) {

    if (params$split_trees) {
        # Split data by mutation type
        data_split <- split(binary_genotypes, binary_genotypes$Mutation_Type)

        for (mutation_type in names(data_split)) {
            # Get mutations from the data
            present_mutations <- get_mutations(
                data, data_split[[mutation_type]]
            )

            # Reassign tree
            reassigned_tree <- reassign_tree(
                tree, present_mutations, params
            )

            # Adjust edge lengths
            reassigned_tree <- adjust_edge_lengths(
                reassigned_tree, params$treemut_pval
            )
            write_and_plot_tree(
                reassigned_tree, params, plot_name = mutation_type
            )
        }
    } else {
        present_mutations <- get_mutations(data, binary_genotypes)
        reassigned_tree <- reassign_tree(
            tree, present_mutations, params
        )
        reassigned_tree <- adjust_edge_lengths(
            reassigned_tree, params$treemut_pval
        )
        write_and_plot_tree(reassigned_tree, params, plot_name = "SNV")
    }
}
