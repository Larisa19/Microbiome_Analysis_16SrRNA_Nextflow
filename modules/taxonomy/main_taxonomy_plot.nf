process TAXONOMY_PLOT {

    tag "Phylum composition"

    publishDir "results/figures", mode: 'copy'

    input:
    path phylum_abundance
    path samplesheet

    output:
    path "phylum_composition.png"

    script:
    """
    Rscript - "$phylum_abundance" "$samplesheet" <<'RSCRIPT'

    args <- commandArgs(trailingOnly = TRUE)

    abundance_file <- args[1]
    samplesheet_file <- args[2]

    library(ggplot2)

    abundance <- read.table(
        abundance_file,
        header = TRUE,
        sep = "\\t",
        row.names = 1,
        check.names = FALSE
    )

    metadata <- read.csv(
        samplesheet_file,
        stringsAsFactors = FALSE
    )

    # Select the 10 most abundant phyla across all samples
    mean_abundance <- rowMeans(abundance)

    top_phyla <- names(
        sort(mean_abundance, decreasing = TRUE)
    )[1:min(10, nrow(abundance))]

    abundance_plot <- abundance

    abundance_plot\$Phylum <- rownames(abundance_plot)

    abundance_plot\$Phylum[
        !abundance_plot\$Phylum %in% top_phyla
    ] <- "Other"

    abundance_plot <- aggregate(
        . ~ Phylum,
        data = abundance_plot,
        FUN = sum
    )

    abundance_long <- reshape(
        abundance_plot,
        varying = colnames(abundance_plot)[-1],
        v.names = "Abundance",
        timevar = "Sample",
        times = colnames(abundance_plot)[-1],
        direction = "long"
    )

    rownames(abundance_long) <- NULL

    abundance_long <- merge(
        abundance_long,
        metadata,
        by.x = "Sample",
        by.y = "sample"
    )

    abundance_long\$treatment <- factor(
        abundance_long\$treatment,
        levels = c("tillage", "cover_crop")
    )

    # Mean abundance by treatment
    treatment_abundance <- aggregate(
        Abundance ~ treatment + Phylum,
        data = abundance_long,
        FUN = mean
    )

    treatment_abundance\$Phylum <- factor(
        treatment_abundance\$Phylum,
        levels = c(top_phyla, "Other")
    )

    p <- ggplot(
        treatment_abundance,
        aes(
            x = treatment,
            y = Abundance,
            fill = Phylum
        )
    ) +
        geom_bar(
            stat = "identity",
            position = "stack"
        ) +
        labs(
            title = "Soil microbial community composition",
            subtitle = "Top 10 phyla — Xynisteri vineyard soil at harvest, non-irrigated",
            x = "Soil management",
            y = "Mean relative abundance (%)",
            fill = "Phylum"
        ) +
        theme_classic()

    ggsave(
        "phylum_composition.png",
        plot = p,
        width = 8,
        height = 6,
        dpi = 300
    )

RSCRIPT
    """

    stub:
    """
    touch phylum_composition.png
    """
}