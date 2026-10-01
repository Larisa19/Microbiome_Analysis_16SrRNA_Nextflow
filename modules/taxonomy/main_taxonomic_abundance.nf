process TAXONOMIC_ABUNDANCE {

    tag "Phylum abundance"

    publishDir "results/taxonomy", mode: 'copy'

    input:
    path asv_table
    path taxonomy

    output:
    path "phylum_abundance.tsv", emit: phylum_abundance

    script:
    """
    Rscript - "$asv_table" "$taxonomy" <<'RSCRIPT'

    args <- commandArgs(trailingOnly = TRUE)

    asv_file <- args[1]
    taxonomy_file <- args[2]

    asv <- read.table(
        asv_file,
        header = TRUE,
        sep = "\\t",
        row.names = 1,
        check.names = FALSE
    )

    tax <- read.table(
        taxonomy_file,
        header = TRUE,
        sep = "\\t",
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    tax <- tax[
        match(rownames(asv), tax\$ASV),
        ,
        drop = FALSE
    ]

    if (!all(rownames(asv) == tax\$ASV)) {
        stop("ASV order does not match between abundance and taxonomy tables.")
    }

    phylum <- tax\$Phylum

    phylum[is.na(phylum) | phylum == ""] <- "Unclassified"

    phylum_abundance <- rowsum(
        as.matrix(asv),
        group = phylum,
        reorder = FALSE
    )

    relative_abundance <- sweep(
        phylum_abundance,
        2,
        colSums(phylum_abundance),
        "/"
    )

    relative_abundance <- relative_abundance * 100

    write.table(
        relative_abundance,
        file = "phylum_abundance.tsv",
        sep = "\\t",
        quote = FALSE,
        col.names = NA
    )

RSCRIPT
    """

    stub:
    """
    touch phylum_abundance.tsv
    """
}