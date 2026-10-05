process DIFFERENTIAL_ABUNDANCE {

    tag "Differential abundance"

    publishDir "results/differential_abundance", mode: 'copy'

    input:
    path asv_table
    path samplesheet
    path taxonomy

    output:
    path "differential_abundance.tsv", emit: results
    path "exploratory_candidates.tsv", emit: candidates
    path "exploratory_candidate_taxa.tsv", emit: candidate_taxa
    path "differential_abundance_summary.txt", emit: summary

    script:

    """
    Rscript - "$asv_table" "$samplesheet" "$taxonomy" <<'RSCRIPT'

    args <- commandArgs(trailingOnly = TRUE)

    asv_file <- args[1]
    metadata_file <- args[2]
    taxonomy_file <- args[3]

    # ------------------------------------------------------------
    # 1. Read input files
    # ------------------------------------------------------------

    asv <- read.table(
        asv_file,
        header = TRUE,
        sep = "\\t",
        row.names = 1,
        check.names = FALSE
    )

    metadata <- read.csv(
        metadata_file,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    taxonomy <- read.table(
        taxonomy_file,
        header = TRUE,
        sep = "\\t",
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    if (!"ASV" %in% colnames(taxonomy)) {
        stop("Taxonomy file must contain an 'ASV' column.")
    }

    # Convert ASV abundance columns explicitly to numeric
    asv[] <- lapply(
        asv,
        function(x) as.numeric(as.character(x))
    )

    # ------------------------------------------------------------
    # 2. Validate metadata
    # ------------------------------------------------------------

    if (!"sample" %in% colnames(metadata)) {
        stop("Metadata must contain a 'sample' column.")
    }

    if (!"treatment" %in% colnames(metadata)) {
        stop("Metadata must contain a 'treatment' column.")
    }

    expected_treatments <- c("tillage", "cover_crop")

    if (!all(expected_treatments %in% metadata\$treatment)) {
        stop(
            "Expected treatment groups 'tillage' and 'cover_crop' were not found."
        )
    }

    # ------------------------------------------------------------
    # 3. Match ASV table and metadata
    # ------------------------------------------------------------

    common_samples <- intersect(
        colnames(asv),
        metadata\$sample
    )

    if (length(common_samples) < 4) {
        stop(
            "Fewer than 4 shared samples were found between ASV table and metadata."
        )
    }

    asv <- asv[, common_samples, drop = FALSE]

    metadata <- metadata[
        match(common_samples, metadata\$sample),
        ,
        drop = FALSE
    ]

    metadata\$treatment <- factor(
        metadata\$treatment,
        levels = c("tillage", "cover_crop")
    )

    # ------------------------------------------------------------
    # 4. Check treatment sample sizes
    # ------------------------------------------------------------

    n_tillage <- sum(metadata\$treatment == "tillage")
    n_cover <- sum(metadata\$treatment == "cover_crop")

    if (n_tillage < 2 || n_cover < 2) {
        stop("At least two samples per treatment are required.")
    }

    # ------------------------------------------------------------
    # 5. Filter low-prevalence ASVs
    # ------------------------------------------------------------

    prevalence_tillage <- rowSums(
        asv[, metadata\$treatment == "tillage", drop = FALSE] > 0
    )

    prevalence_cover <- rowSums(
        asv[, metadata\$treatment == "cover_crop", drop = FALSE] > 0
    )

    keep <- (
        prevalence_tillage >= 2 |
        prevalence_cover >= 2
    )

    asv_filtered <- asv[keep, , drop = FALSE]

    # ------------------------------------------------------------
    # 6. Convert counts to relative abundance
    # ------------------------------------------------------------

    sample_totals <- colSums(asv_filtered)

    if (any(sample_totals == 0)) {
        stop(
            "At least one sample has zero reads after prevalence filtering."
        )
    }

    relative_abundance <- sweep(
        asv_filtered,
        2,
        sample_totals,
        "/"
    )

    # ------------------------------------------------------------
    # 7. Differential abundance testing
    # ------------------------------------------------------------

    results <- data.frame(
        ASV = rownames(relative_abundance),
        prevalence_tillage = prevalence_tillage[keep],
        prevalence_cover_crop = prevalence_cover[keep],
        mean_tillage = NA_real_,
        mean_cover_crop = NA_real_,
        median_tillage = NA_real_,
        median_cover_crop = NA_real_,
        log2_fold_change = NA_real_,
        p_value = NA_real_,
        stringsAsFactors = FALSE
    )

    tillage_samples <- metadata\$treatment == "tillage"
    cover_samples <- metadata\$treatment == "cover_crop"

    for (i in seq_len(nrow(relative_abundance))) {

        x <- as.numeric(
            relative_abundance[i, tillage_samples]
        )

        y <- as.numeric(
            relative_abundance[i, cover_samples]
        )

        results\$mean_tillage[i] <- mean(x)
        results\$mean_cover_crop[i] <- mean(y)

        results\$median_tillage[i] <- median(x)
        results\$median_cover_crop[i] <- median(y)

        pseudocount <- 1e-6

        results\$log2_fold_change[i] <- log2(
            (mean(y) + pseudocount) /
            (mean(x) + pseudocount)
        )

        test <- suppressWarnings(
            wilcox.test(
                x,
                y,
                exact = FALSE
            )
        )

        results\$p_value[i] <- test\$p.value
    }

    # ------------------------------------------------------------
    # 8. Multiple-testing correction
    # ------------------------------------------------------------

    results\$FDR <- p.adjust(
        results\$p_value,
        method = "BH"
    )

    # ------------------------------------------------------------
    # 9. Classify direction of change
    # ------------------------------------------------------------

    results\$direction <- "No clear difference"

    results\$direction[
        results\$log2_fold_change > 0
    ] <- "Higher in cover_crop"

    results\$direction[
        results\$log2_fold_change < 0
    ] <- "Higher in tillage"

    # ------------------------------------------------------------
    # 10. Flag statistically significant ASVs
    # ------------------------------------------------------------

    results\$significant <- (
        results\$FDR < 0.05 &
        abs(results\$log2_fold_change) >= 1
    )

    # ------------------------------------------------------------
    # 10b. Flag exploratory candidates
    # ------------------------------------------------------------

    results\$exploratory_candidate <- (
        abs(results\$log2_fold_change) >= 2 &
        (
            results\$prevalence_tillage >= 2 |
            results\$prevalence_cover_crop >= 2
        )
    )

    # ------------------------------------------------------------
    # 11. Order results
    # ------------------------------------------------------------

    results <- results[
        order(
            results\$FDR,
            -abs(results\$log2_fold_change)
        ),
    ]

    # ------------------------------------------------------------
    # 12. Write differential abundance results
    # ------------------------------------------------------------

    write.table(
        results,
        file = "differential_abundance.tsv",
        sep = "\\t",
        quote = FALSE,
        row.names = FALSE
    )

    # ------------------------------------------------------------
    # 12b. Write exploratory candidates
    # ------------------------------------------------------------

    exploratory_candidates <- results[
        results\$exploratory_candidate,
        ,
        drop = FALSE
    ]

    write.table(
        exploratory_candidates,
        file = "exploratory_candidates.tsv",
        sep = "\\t",
        quote = FALSE,
        row.names = FALSE
    )

    # ------------------------------------------------------------
    # 12c. Add taxonomy to exploratory candidates
    # ------------------------------------------------------------

    exploratory_candidate_taxa <- merge(
        exploratory_candidates,
        taxonomy,
        by = "ASV",
        all.x = TRUE
    )

    taxonomy_columns <- c(
        "ASV",
        "Kingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus"
    )

    existing_taxonomy_columns <- taxonomy_columns[
        taxonomy_columns %in% colnames(exploratory_candidate_taxa)
    ]

    exploratory_candidate_taxa <- exploratory_candidate_taxa[
        ,
        c(
            existing_taxonomy_columns,
            setdiff(
                colnames(exploratory_candidate_taxa),
                existing_taxonomy_columns
            )
        ),
        drop = FALSE
    ]

    write.table(
        exploratory_candidate_taxa,
        file = "exploratory_candidate_taxa.tsv",
        sep = "\\t",
        quote = FALSE,
        row.names = FALSE
    )

    # ------------------------------------------------------------
    # 13. Summary
    # ------------------------------------------------------------

    n_total <- nrow(asv)
    n_filtered <- nrow(asv_filtered)
    n_candidates <- sum(results\$significant)

    n_cover_candidates <- sum(
        results\$significant &
        results\$direction == "Higher in cover_crop"
    )

    n_tillage_candidates <- sum(
        results\$significant &
        results\$direction == "Higher in tillage"
    )

    n_exploratory <- sum(
        results\$exploratory_candidate
    )

    summary_lines <- c(
        "Differential abundance analysis",
        "================================",
        "",
        paste("Samples analysed:", length(common_samples)),
        paste("Tillage samples:", n_tillage),
        paste("Cover-crop samples:", n_cover),
        "",
        paste("Total ASVs:", n_total),
        paste(
            "ASVs retained after prevalence filtering:",
            n_filtered
        ),
        "",
        "Statistical test: Wilcoxon rank-sum test",
        "Multiple-testing correction: Benjamini-Hochberg FDR",
        "Significance threshold: FDR < 0.05",
        "Significant effect-size threshold: |log2 fold change| >= 1",
        "Exploratory candidate threshold: |log2 fold change| >= 2",
        "Exploratory prevalence threshold: >= 2 samples in at least one treatment",
        "",
        paste(
            "Statistically significant ASVs:",
            n_candidates
        ),
        paste(
            "Higher in cover_crop:",
            n_cover_candidates
        ),
        paste(
            "Higher in tillage:",
            n_tillage_candidates
        ),
        paste(
            "Exploratory candidate ASVs:",
            n_exploratory
        ),
        "",
        "Interpretation:",
        "Results are exploratory because only four samples were available per treatment."
    )

    writeLines(
        summary_lines,
        con = "differential_abundance_summary.txt"
    )

    RSCRIPT
    """
}