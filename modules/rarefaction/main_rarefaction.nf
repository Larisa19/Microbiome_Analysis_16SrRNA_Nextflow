process RAREFACTION {

    tag "Rarefaction"

    publishDir "results/rarefaction", mode: 'copy'

    input:
    path asv_table

    output:
    path "asv_table_rarefied.tsv", emit: rarefied_table
    path "rarefaction_info.txt", emit: rarefaction_info

    script:
    """
    Rscript - "$asv_table" <<'RSCRIPT'

    args <- commandArgs(trailingOnly = TRUE)

    asv_file <- args[1]

    # Read ASV table
    asv <- read.table(
        asv_file,
        header = TRUE,
        sep = "\\t",
        row.names = 1,
        check.names = FALSE
    )

    # Calculate sequencing depth per sample
    sample_depths <- colSums(asv)

    # Use the minimum sequencing depth across samples
    rarefaction_depth <- min(sample_depths)

    # Rarefy each sample without replacement
    set.seed(123)

    rarefy_sample <- function(counts, depth) {

        if (sum(counts) < depth) {
            stop("Sample has fewer reads than the rarefaction depth.")
        }

        reads <- rep(seq_along(counts), counts)

        selected <- sample(
            reads,
            size = depth,
            replace = FALSE
        )

        rarefied <- tabulate(
            selected,
            nbins = length(counts)
        )

        return(rarefied)
    }

    rarefied <- as.data.frame(
        lapply(asv, rarefy_sample, depth = rarefaction_depth)
    )

    rownames(rarefied) <- rownames(asv)

    # Save rarefied ASV table
    write.table(
        rarefied,
        file = "asv_table_rarefied.tsv",
        sep = "\\t",
        quote = FALSE,
        col.names = NA
    )

    # Save information for reproducibility
    info <- data.frame(
        sample = names(sample_depths),
        original_reads = as.integer(sample_depths),
        rarefaction_depth = rarefaction_depth
    )

    write.table(
        info,
        file = "rarefaction_info.txt",
        sep = "\\t",
        quote = FALSE,
        row.names = FALSE
    )

    RSCRIPT
    """

    stub:
    """
    touch asv_table_rarefied.tsv
    touch rarefaction_info.txt
    """
}