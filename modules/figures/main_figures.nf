process ALPHA_DIVERSITY_PLOT {

    tag "alpha_diversity"

    publishDir "results/figures", mode: 'copy'

    input:
    path alpha_table

    output:
    path "alpha_diversity_shannon.png"

    script:
    """
    Rscript - <<'EOF'

    library(ggplot2)

    alpha <- read.table(
        "$alpha_table",
        header = TRUE,
        sep = "\\t",
        stringsAsFactors = FALSE
    )

    alpha\$treatment <- factor(
        alpha\$treatment,
        levels = c("tillage", "cover_crop")
    )

    p <- ggplot(
        alpha,
        aes(x = treatment, y = Shannon)
    ) +
        geom_boxplot(
            width = 0.5,
            outlier.shape = NA
        ) +
        geom_jitter(
            width = 0.08,
            size = 3
        ) +
        labs(
            title = "Alpha diversity of soil microbial communities",
            subtitle = "Xynisteri vineyard soil at harvest — non-irrigated",
            x = "Soil management",
            y = "Shannon diversity"
        ) +
        theme_classic(base_size = 13)

    ggsave(
        "alpha_diversity_shannon.png",
        plot = p,
        width = 7,
        height = 5,
        dpi = 300
    )

    EOF
    """

    stub:
    """
    touch alpha_diversity_shannon.png
    """
}
