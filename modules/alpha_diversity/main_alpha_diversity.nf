process ALPHA_DIVERSITY {

    tag "Alpha diversity"

    publishDir "results/alpha_diversity", mode: 'copy'

    input:
    path asv_table
    path samplesheet

    output:
    path "alpha_diversity.tsv", emit: alpha_diversity

    script:

    """
    Rscript - "$asv_table" "$samplesheet" <<'RSCRIPT'

    args <- commandArgs(trailingOnly = TRUE)

    asv_file <- args[1]
    metadata_file <- args[2]

    # Read ASV table
    asv <- read.table(
        asv_file,
        header = TRUE,
        row.names = 1,
        sep = "\\t",
        check.names = FALSE
    )

    # Read metadata
    metadata <- read.csv(
        metadata_file,
        stringsAsFactors = FALSE
    )

    # Validate sample IDs
    missing_metadata <- setdiff(colnames(asv), metadata\$sample)

    if (length(missing_metadata) > 0) {
        stop(
            paste(
                "Samples missing from metadata:",
                paste(missing_metadata, collapse = ", ")
            )
        )
    }

    # Calculate observed ASVs
    observed_asvs <- colSums(asv > 0)

    # Calculate Shannon diversity
    shannon <- apply(asv, 2, function(x) {
        if (sum(x) == 0) {
            return(0)
        }

        p <- x[x > 0] / sum(x)

        -sum(p * log(p))
    })

    # Build results
    results <- data.frame(
        sample = colnames(asv),
        observed_ASVs = as.numeric(observed_asvs),
        Shannon = as.numeric(shannon)
    )

    # Add metadata
    results <- merge(
        metadata,
        results,
        by = "sample",
        all.y = TRUE,
        sort = FALSE
    )

    # Restore ASV table order
    results <- results[match(colnames(asv), results\$sample), ]

    # Write results
    write.table(
        results,
        file = "alpha_diversity.tsv",
        sep = "\\t",
        row.names = FALSE,
        quote = FALSE
    )

    RSCRIPT
    """

    stub:
    """
    touch alpha_diversity.tsv
    """
}
