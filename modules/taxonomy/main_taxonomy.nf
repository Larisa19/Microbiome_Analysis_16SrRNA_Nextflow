process TAXONOMY {

    tag "ASV taxonomy"

    publishDir "results/taxonomy", mode: 'copy'

    input:
    path asv_table
    path reference

    output:
    path "taxonomy.tsv", emit: taxonomy

    script:
    """
    Rscript - "$asv_table" "$reference" <<'RSCRIPT'

    args <- commandArgs(trailingOnly = TRUE)

    asv_file <- args[1]
    reference_file <- args[2]

    if (!requireNamespace("dada2", quietly = TRUE)) {
        stop("R package 'dada2' is required.")
    }

    asv <- read.table(
        asv_file,
        header = TRUE,
        sep = "\\t",
        row.names = 1,
        check.names = FALSE
    )

    seqs <- rownames(asv)

    taxonomy <- dada2::assignTaxonomy(
        seqs,
        reference_file,
        multithread = TRUE
    )

    taxonomy <- as.data.frame(taxonomy)

    taxonomy\$ASV <- seqs

    taxonomy <- taxonomy[
        ,
        c(
            "ASV",
            setdiff(
                colnames(taxonomy),
                "ASV"
            )
        )
    ]

    write.table(
        taxonomy,
        file = "taxonomy.tsv",
        sep = "\\t",
        quote = FALSE,
        row.names = FALSE
    )

RSCRIPT
    """

    stub:
    """
    touch taxonomy.tsv
    """
}