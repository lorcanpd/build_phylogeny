library(VariantAnnotation)
library(dplyr)
library(data.table)
library(ggplot2)


# function to read Caveman file
# output df with Chr_Pos_Ref_Alt, NV, NR, VAF=NV/NR, Gene
# filter for PASS, ASMD, CLMP
caveman2 <- function(sample, path, filter = TRUE, ASMD_value = 140, CLMP_value = 0) {

    data <- readVcf(paste0(path, "/", sample, ".caveman_c.annot.vcf.gz"))
    # Number of reads
    reads <- cbind(
        (geno(data)$FAZ)[, 2] + (geno(data)$RAZ)[, 2],
        (geno(data)$FCZ)[, 2] + (geno(data)$RCZ)[, 2],
        (geno(data)$FGZ)[, 2] + (geno(data)$RGZ)[, 2],
        (geno(data)$FTZ)[, 2] + (geno(data)$RTZ)[, 2]
    )
    Var <- data.frame(
        Sample = sample,
        Chr = as.character(seqnames(rowRanges(data))),
        Pos = start(ranges(rowRanges(data))),
        Ref = as.character(ref(data))
    )
    Alt_tmp <- CharacterList(alt(data))
    # check all Variants have an alt allele
    Delete <- which(
        sum(alt(data) %in% c("A", "C", "T", "G")) != 1
    )
    Alt_tmp[Delete] <- "N"
    Var$Alt <- as.character(unlist(Alt_tmp))
    Var$NR <- rowSums(reads) #Total number of reads
    Var$NV <- NA
    Var$Gene <- info(data)$VD
    Var$Depth <- info(data)$DP
    Var$Gene <- gsub("\\|.*", "", Var$Gene)
    Var$Impact <- info(data)$VC
    Var$tmp <- info(data)$VW
    Var$Snp <- info(data)$SNP
    Var$AAchange <- gsub(".*p\\.", "", Var$tmp)
    Var$AAchange <- gsub("\\|.*", "", Var$AAchange)
    colnames(reads) <- c("A", "C", "G", "T")

    for (k in c("A", "C", "G", "T")) {
        # number of reads with the alternate base
        Var$NV[Var$Alt == k] <- reads[Var$Alt == k, k]
    }

    if (filter) {
        select <- rowRanges(data)$FILTER == "PASS"
        Var$Flag <- info(data)$ASMD >= ASMD_value & info(data)$CLPM == CLMP_value
        Var <- Var[select, ]
    }

    Var$ID <- paste(Var$Chr, Var$Pos, Var$Ref, Var$Alt, sep = "_")
    Var$VAF <- Var$NV / Var$NR
    return(Var)
}

# Function to count unique and shared variants
count_unique_shared <- function(data) {
    # Get list of unique variants per sample
    unique_variants <- lapply(
        split(data, data$Sample),
        function(x) {
            unique(x$ID)
        }
    )
    # Get list of shared variants per sample
    shared_variants <- lapply(
        split(data, data$Sample),
        function(x) {
            unique(x$ID[duplicated(x$ID)])
        }
    )
    # Get list of unique variants per sample
    unique_counts <- lapply(unique_variants, length)
    # Get list of shared variants per sample
    shared_counts <- lapply(shared_variants, length)
    # Combine into data frame
    counts <- data.frame(
        Sample = names(unique_counts),
        Unique = unlist(unique_counts),
        Shared = unlist(shared_counts),
        stringsAsFactors = FALSE
    )
    return(counts)
}

# path to samples directory
path <- "/nfs/cancer_ref01/nst_links/live/3187"

# get all directories in path
dirs <- list.dirs(path, recursive = FALSE)

# Create empty data frame to store all samples, using column names from above
all_samples <- data.frame(
    Sample = character(),
    Chr = character(),
    Pos = numeric(),
    Ref = character(),
    Alt = character(),
    NV = numeric(),
    NR = numeric(),
    Gene = character(),
    Depth = numeric(),
    Impact = character(),
    Snp = character(),
    AAchange = character(),
    Flag = logical(),
    ID = character(),
    VAF = numeric(),
    stringsAsFactors = FALSE
)

for (dir in dirs) {
    # get sample name from directory name
    sample <- basename(dir)
    # read in file
    data <- caveman2(sample = sample, path = dir)
    # add to all samples
    all_samples <- rbind(all_samples, data)
}

# Count unique and shared variants
counts <- count_unique_shared(all_samples)

# Do barplot of unique and shared variants per sample, with the unique variants
# stacked on top of the shared variants for each sample
counts_long <- counts %>%
    gather(key = "Type", value = "Count", Unique, Shared) %>%
    mutate(Type = factor(Type, levels = c("Unique", "Shared")))

ggplot(counts_long, aes(x = Sample, y = Count, fill = Type)) +
    geom_bar(stat = "identity", position = "stack") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Sample", y = "Number of variants", fill = "Type")


