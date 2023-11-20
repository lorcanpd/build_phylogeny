
require(c("Rsamtools", "GenomicRanges", "dplyr"))

plot_spectrum <- function(
    bed, save, add_to_title="",
    genomeFile = "/nfs/cancer_ref01/Homo_sapiens/37/genome.fa"
){
    mutations <- as.data.frame(bed)
    colnames(mutations) <- c("chr", "pos", "ref", "mut")
    mutations$pos <- as.numeric(mutations$pos)
    mutations <- mutations[
        (mutations$ref %in% c("A", "C", "G", "T")) &
            (mutations$mut %in% c("A", "C", "G", "T")) &
            (mutations$chr %in% c(1:22, "X", "Y")),
    ]
    # Extract genomic ranges using given ranges determined by GRanges and
    # IRanges
    mutations$trinuc_ref <- as.vector(
        scanFa(
            genomeFile,
            GRanges(
                mutations$chr,
                IRanges(as.numeric(mutations$pos) - 1,
                        as.numeric(mutations$pos) + 1)
            )
        )
    )
    # 2. Annotating the mutation from the pyrimidine base
    ntcomp <- c(T="A", G="C", C="G", A="T")
    mutations <- mutations %>%
        mutate(
            sub = if_else(
                ref %in% c("A", "G"),
                paste(ntcomp[ref], ntcomp[mut], sep=">"),
                paste(ref, mut, sep=">")
            ),
            trinuc_ref_py = if_else(
                ref %in% c("A", "G"),
                paste(
                    ntcomp[rev(strsplit(trinuc_ref, split="")[[1]])],
                    collapse=""
                ),
                trinuc_ref
            )
        )
    # REPLACED WITH DPLYR MUTATE
    # mutations$sub <- paste(mutations$ref, mutations$mut, sep=">")
    # mutations$trinuc_ref_py <- mutations$trinuc_ref
    # for (j in 1:nrow(mutations)) {
    #     if (mutations$ref[j] %in% c("A","G")) { # Purine base
    #         mutations$sub[j] <- paste(ntcomp[mutations$ref[j]], ntcomp[mutations$mut[j]], sep=">")
    #         mutations$trinuc_ref_py[j] <- paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j], split="")[[1]])], collapse="")
    #     }
    # }


    # 3. Counting subs
    freqs <- table(
        paste(
            mutations$sub,
            paste(
                substr(mutations$trinuc_ref_py, 1, 1),
                substr(mutations$trinuc_ref_py, 3, 3),
                sep="-"
            ),
            sep=","
        )
    )
    sub_vec <- c("C>A","C>G","C>T","T>A","T>C","T>G")
    ctx_vec <- paste(
        rep(c("A","C","G","T"), each=4),
        rep(c("A","C","G","T"), times=4),
        sep="-"
    )
    full_vec <- paste(rep(sub_vec,each=16), rep(ctx_vec,times=6),sep=",")
    freqs_full <- freqs[full_vec]
    freqs_full[is.na(freqs_full)] <- 0
    names(freqs_full) <- full_vec

    xstr <- paste(
        substr(full_vec,5,5),
        substr(full_vec,1,1),
        substr(full_vec,7,7),
        sep=""
    )

    if(!is.null(save)) {
        pdf(save, width = 12, height = 4)
    }
    if(is.null(save)) {
        dev.new(width = 12, height = 4)
    }
    colvec <- rep(
        c("dodgerblue", "black", "red", "grey70", "olivedrab3", "plum2"),
        each=16
    )
    y <- freqs_full
    maxy <- max(y)
    h <- barplot(y, las=2, col=colvec, border=NA, ylim=c(0,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, ylab="Number mutations", main=paste0("Number of mutations: ",sum(freqs_full), add_to_title))
    for (j in 1:length(sub_vec)) {
        xpos <- h[c((j-1)*16+1,j*16)]
        rect(
            xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA,
            col=colvec[j*16]
        )
        text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
    }
    if(!is.null(save)) {
        dev.off()
    }
}


library(dplyr)

exact.binomial <- function(gender, NV, NR, cutoff=-5, qval_return=F) {
    # Function to filter out germline variants based on unmatched
    # variant calls of multiple samples from same individual (aggregate coverage
    # ideally >150 or so, but will work with less). NV is matrix of reads
    # supporting variants and NR the matrix with total depth (samples as
    # columns, mutations rows, with rownames as chr_pos_ref_alt or equivalent).
    # Returns a logical vector, TRUE if mutation is likely to be germline.

    # Convert NV and NR to data frames for dplyr processing
    NV_df <- as.data.frame(NV)
    NR_df <- as.data.frame(NR)

    # Combine NV and NR with their row names
    combined_df <- NV_df %>%
        mutate(chr_pos = rownames(NV_df)) %>%
        rowwise() %>%
        mutate(
            NV_sum = sum(c_across(-chr_pos)),
            NR_sum = sum(NR_df[chr_pos, ])
        )

    XY_chromosomal <- grepl("X|Y", combined_df$chr_pos)
    autosomal <- !XY_chromosomal

    # Perform binomial test
    combined_df <- combined_df %>%
        mutate(
            pval = case_when(
                gender == "female" & NR_sum > 0 ~
                    binom.test(
                        x = NV_sum, n = NR_sum, p = 0.5, alt = 'less'
                    )$p.value,
                gender == "male" & autosomal & NR_sum > 0 ~
                    binom.test(
                        x = NV_sum, n = NR_sum, p = 0.5, alt = 'less'
                    )$p.value,
                gender == "male" & XY_chromosomal & NR_sum > 0 ~
                    binom.test(
                        x = NV_sum, n = NR_sum, p = 0.95, alt = 'less'
                    )$p.value,
                TRUE ~ 1
            )
        )

    # Adjust p-values and determine germline status
    qval <- p.adjust(combined_df$pval, method = "BH")

    if (qval_return) {
        return(qval)
    } else {
        germline <- log10(qval) > cutoff
        return(germline)
    }
}

# REPLACED WITH DPLYR VERSION
# exact.binomial <- function(gender, NV, NR, cutoff=-5, qval_return=F) {
#     # Function to filter out germline variants based on unmatched
#     # variant calls of multiple samples from same individual (aggregate coverage
#     # ideally >150 or so, but will work with less). NV is matrix of reads
#     # supporting variants and NR the matrix with total depth (samples as
#     # columns, mutations rows, with rownames as chr_pos_ref_alt or equivalent).
#     # Returns a logical vector, TRUE if mutation is likely to be germline.
#
#     XY_chromosomal <- grepl("X|Y", rownames(NR))
#     autosomal <- !XY_chromosomal
#
#     if (gender == "female") {
#         NV_vec <- rowSums(NV)
#         NR_vec <- rowSums(NR)
#         pval <- rep(1, length(NV_vec))
#         for (n in 1:length(NV_vec)) {
#             if (NR_vec[n]>0){
#                 pval[n] <- binom.test(
#                     x=NV_vec[n],
#                     n=NR_vec[n],
#                     p=0.5,
#                     alt='less'
#                 )$p.value
#             }
#         }
#     }
#     # For male, split test in autosomal and XY chromosomal part
#     if (gender == "male") {
#         pval <- rep(1,nrow(NV))
#         NV_vec <- rowSums(NV)[autosomal]
#         NR_vec <- rowSums(NR)[autosomal]
#         pval_auto <- rep(1, sum(autosomal))
#         pval_XY <- rep(1, sum(XY_chromosomal))
#
#         for (n in 1:sum(autosomal)){
#             if (NR_vec[n]>0) {
#                 pval_auto[n] <- binom.test(
#                     x=NV_vec[n],
#                     n=NR_vec[n],
#                     p=0.5,
#                     alt='less'
#                 )$p.value
#             }
#         }
#         NV_vec <- rowSums(NV)[XY_chromosomal]
#         NR_vec <- rowSums(NR)[XY_chromosomal]
#         for (n in 1:sum(XY_chromosomal)) {
#             if (NR_vec[n] > 0) {
#                 pval_XY[n] <- binom.test(
#                     x=NV_vec[n],
#                     n=NR_vec[n],
#                     p=0.95,
#                     alt='less'
#                 )$p.value
#             }
#         }
#         pval[autosomal] <- pval_auto
#         pval[XY_chromosomal] <- pval_XY
#     }
#     qval <- p.adjust(pval, method="BH")
#     if (qval_return) {
#         return(qval)
#     } else {
#         germline <- log10(qval) > cutoff
#         return(germline)
#     }
# }