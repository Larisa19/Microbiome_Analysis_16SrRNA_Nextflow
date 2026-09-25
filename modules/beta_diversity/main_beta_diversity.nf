process BETA_DIVERSITY {

    tag "Beta diversity"

    publishDir "results/beta_diversity", mode: 'copy'

    input:
    path rarefied_table
    path samplesheet

    output:
    path "bray_curtis_distance.tsv", emit: bray_curtis
    path "pcoa_coordinates.tsv", emit: pcoa_coordinates
    path "permanova_results.tsv", emit: permanova
    path "beta_diversity_summary.txt", emit: summary
    path "pcoa_treatment.png", emit: pcoa_plot

    script:
    """
    Rscript - "$rarefied_table" "$samplesheet" <<'RSCRIPT'

    args <- commandArgs(trailingOnly = TRUE)

    rarefied_file <- args[1]
    samplesheet_file <- args[2]

    # Required packages
    if (!requireNamespace("vegan", quietly = TRUE)) {
        stop("R package 'vegan' is required.")
    }

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("R package 'ggplot2' is required.")
    }

    # Read rarefied ASV table
    asv <- read.table(
        rarefied_file,
        header = TRUE,
        sep = "\\t",
        row.names = 1,
        check.names = FALSE
    )

    # Read metadata
    metadata <- read.csv(
        samplesheet_file,
        stringsAsFactors = FALSE
    )

    # Check sample IDs
    missing_metadata <- setdiff(
        colnames(asv),
        metadata\$sample
    )

    if (length(missing_metadata) > 0) {
        stop(
            paste(
                "Samples missing from metadata:",
                paste(missing_metadata, collapse = ", ")
            )
        )
    }

    # Keep metadata in ASV-table order
    metadata <- metadata[
        match(
            colnames(asv),
            metadata\$sample
        ),
        ,
        drop = FALSE
    ]

    # Transpose ASV table:
    # rows = samples, columns = ASVs
    community <- t(asv)

    # Bray-Curtis dissimilarity
    bray <- vegan::vegdist(
        community,
        method = "bray"
    )

    # Save Bray-Curtis matrix
    bray_matrix <- as.matrix(bray)

    write.table(
        bray_matrix,
        file = "bray_curtis_distance.tsv",
        sep = "\\t",
        quote = FALSE,
        col.names = NA
    )

    # PCoA using classical multidimensional scaling
    pcoa <- cmdscale(
        bray,
        k = 2,
        eig = TRUE
    )

    pcoa_coordinates <- data.frame(
        sample = rownames(pcoa\$points),
        PCoA1 = pcoa\$points[, 1],
        PCoA2 = pcoa\$points[, 2]
    )

    # Add metadata
    pcoa_coordinates <- merge(
        pcoa_coordinates,
        metadata,
        by = "sample",
        sort = FALSE
    )

    # Restore ASV-table sample order
    pcoa_coordinates <- pcoa_coordinates[
        match(
            colnames(asv),
            pcoa_coordinates\$sample
        ),
        ,
        drop = FALSE
    ]

    # Save PCoA coordinates
    write.table(
        pcoa_coordinates,
        file = "pcoa_coordinates.tsv",
        sep = "\\t",
        quote = FALSE,
        row.names = FALSE
    )

    # PERMANOVA
    permanova <- vegan::adonis2(
        bray ~ treatment,
        data = metadata,
        permutations = 999
    )

    # Save PERMANOVA results
    permanova_table <- as.data.frame(permanova)

    permanova_table\$term <- rownames(permanova_table)

    permanova_table <- permanova_table[
        ,
        c(
            "term",
            setdiff(
                colnames(permanova_table),
                "term"
            )
        )
    ]

    write.table(
        permanova_table,
        file = "permanova_results.tsv",
        sep = "\\t",
        quote = FALSE,
        row.names = FALSE
    )

    # Variance explained
    variance_explained <- pcoa\$eig / sum(pcoa\$eig)

    # PCoA plot
    pcoa_plot <- ggplot2::ggplot(
        pcoa_coordinates,
        ggplot2::aes(
            x = PCoA1,
            y = PCoA2,
            color = treatment,
            label = sample
        )
    ) +
        ggplot2::geom_point(
            size = 4
        ) +
        ggplot2::geom_text(
            vjust = -0.8,
            size = 3,
            show.legend = FALSE
        ) +
        ggplot2::labs(
            title = "PCoA of soil microbial communities",
            subtitle = "Xynisteri vineyard soil at harvest — non-irrigated",
            x = sprintf(
                "PCoA1 (%.2f%%)",
                variance_explained[1] * 100
            ),
            y = sprintf(
                "PCoA2 (%.2f%%)",
                variance_explained[2] * 100
            ),
            color = "Treatment"
        ) +
        ggplot2::theme_classic()

    # Save PCoA figure
    ggplot2::ggsave(
        "pcoa_treatment.png",
        plot = pcoa_plot,
        width = 8,
        height = 6,
        dpi = 300
    )

    # Summary
    sink("beta_diversity_summary.txt")

    cat("BETA DIVERSITY ANALYSIS\\n")
    cat("=======================\\n\\n")

    cat("Distance metric: Bray-Curtis\\n")
    cat("Ordination: PCoA\\n")
    cat("Statistical test: PERMANOVA\\n")
    cat("Permutations: 999\\n\\n")

    cat("PCoA variance explained:\\n")

    cat(
        sprintf(
            "PCoA1: %.2f%%\\n",
            variance_explained[1] * 100
        )
    )

    cat(
        sprintf(
            "PCoA2: %.2f%%\\n\\n",
            variance_explained[2] * 100
        )
    )

    cat("PERMANOVA:\\n")
    print(permanova)

    sink()

    RSCRIPT
    """

    stub:
    """
    touch bray_curtis_distance.tsv
    touch pcoa_coordinates.tsv
    touch permanova_results.tsv
    touch beta_diversity_summary.txt
    touch pcoa_treatment.png
    """
}