process CANDIDATE_TAXA {

    tag "Exploratory candidate taxa"

    publishDir "results/candidate_taxa", mode: 'copy'

    input:
    path asv_table
    path exploratory_candidates
    path taxonomy
    path samplesheet
    path beta_summary

    output:
    path "candidate_taxa.tsv", emit: candidate_taxa
    path "candidate_taxa_abundance.tsv", emit: abundance
    path "candidate_taxa_heatmap.png", emit: heatmap
    path "candidate_taxa_summary.txt", emit: summary

    script:

    """
    Rscript - "$asv_table" "$exploratory_candidates" "$taxonomy" "$samplesheet" "$beta_summary" <<'RSCRIPT'

    args <- commandArgs(trailingOnly = TRUE)

    asv_file <- args[1]
    candidates_file <- args[2]
    taxonomy_file <- args[3]
    metadata_file <- args[4]
    beta_file <- args[5]

    suppressPackageStartupMessages({
        library(ggplot2)
    })

    # ------------------------------------------------------------
    # Read input files
    # ------------------------------------------------------------

    asv <- read.table(
        asv_file,
        header = TRUE,
        sep = "\\t",
        row.names = 1,
        check.names = FALSE
    )

    candidates <- read.table(
        candidates_file,
        header = TRUE,
        sep = "\\t",
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

    metadata <- read.csv(
        metadata_file,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    # ------------------------------------------------------------
    # Validate input columns
    # ------------------------------------------------------------

    required_candidate_columns <- c(
        "ASV",
        "log2_fold_change",
        "p_value",
        "FDR",
        "direction"
    )

    missing_candidate_columns <- setdiff(
        required_candidate_columns,
        colnames(candidates)
    )

    if (length(missing_candidate_columns) > 0) {
        stop(
            paste(
                "Missing columns in exploratory_candidates.tsv:",
                paste(missing_candidate_columns, collapse = ", ")
            )
        )
    }

    if (!"ASV" %in% colnames(taxonomy)) {
        stop("Taxonomy file must contain an 'ASV' column.")
    }

    if (!"sample" %in% colnames(metadata)) {
        stop("Metadata must contain a 'sample' column.")
    }

    if (!"treatment" %in% colnames(metadata)) {
        stop("Metadata must contain a 'treatment' column.")
    }

    # ------------------------------------------------------------
    # Select exploratory candidates
    # ------------------------------------------------------------

    candidates <- candidates[
        abs(candidates[["log2_fold_change"]]) >= 2,
        ,
        drop = FALSE
    ]

    if (nrow(candidates) == 0) {
        stop(
            "No exploratory candidate ASVs met the |log2 fold change| >= 2 threshold."
        )
    }

    # ------------------------------------------------------------
    # Add taxonomy
    # ------------------------------------------------------------

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
        taxonomy_columns %in% colnames(taxonomy)
    ]

    candidate_taxa <- merge(
        candidates,
        taxonomy[
            ,
            existing_taxonomy_columns,
            drop = FALSE
        ],
        by = "ASV",
        all.x = TRUE
    )

    # ------------------------------------------------------------
    # Create taxonomic label
    # Prefer Genus, then Family, Order, Phylum, then ASV
    # ------------------------------------------------------------

    candidate_taxa[["taxon_label"]] <- candidate_taxa[["ASV"]]

    if ("Genus" %in% colnames(candidate_taxa)) {

        use_genus <- (
            !is.na(candidate_taxa[["Genus"]]) &
            candidate_taxa[["Genus"]] != "" &
            candidate_taxa[["Genus"]] != "NA"
        )

        candidate_taxa[["taxon_label"]][use_genus] <-
            candidate_taxa[["Genus"]][use_genus]
    }

    if ("Family" %in% colnames(candidate_taxa)) {

        use_family <- (
            candidate_taxa[["taxon_label"]] == candidate_taxa[["ASV"]] &
            !is.na(candidate_taxa[["Family"]]) &
            candidate_taxa[["Family"]] != "" &
            candidate_taxa[["Family"]] != "NA"
        )

        candidate_taxa[["taxon_label"]][use_family] <-
            candidate_taxa[["Family"]][use_family]
    }

    if ("Order" %in% colnames(candidate_taxa)) {

        use_order <- (
            candidate_taxa[["taxon_label"]] == candidate_taxa[["ASV"]] &
            !is.na(candidate_taxa[["Order"]]) &
            candidate_taxa[["Order"]] != "" &
            candidate_taxa[["Order"]] != "NA"
        )

        candidate_taxa[["taxon_label"]][use_order] <-
            candidate_taxa[["Order"]][use_order]
    }

    if ("Phylum" %in% colnames(candidate_taxa)) {

        use_phylum <- (
            candidate_taxa[["taxon_label"]] == candidate_taxa[["ASV"]] &
            !is.na(candidate_taxa[["Phylum"]]) &
            candidate_taxa[["Phylum"]] != "" &
            candidate_taxa[["Phylum"]] != "NA"
        )

        candidate_taxa[["taxon_label"]][use_phylum] <-
            candidate_taxa[["Phylum"]][use_phylum]
    }

    # ------------------------------------------------------------
    # Rank candidates by effect size
    # ------------------------------------------------------------

    candidate_taxa <- candidate_taxa[
        order(
            -abs(candidate_taxa[["log2_fold_change"]]),
            candidate_taxa[["FDR"]]
        ),
        ,
        drop = FALSE
    ]

    # Keep the ten strongest exploratory candidates
    max_candidates <- min(10, nrow(candidate_taxa))

    candidate_taxa <- candidate_taxa[
        seq_len(max_candidates),
        ,
        drop = FALSE
    ]

    # ------------------------------------------------------------
    # Match candidates to ASV table
    # ------------------------------------------------------------

    common_asvs <- intersect(
        candidate_taxa[["ASV"]],
        rownames(asv)
    )

    if (length(common_asvs) == 0) {
        stop("None of the candidate ASVs were found in the ASV table.")
    }

    asv <- asv[
        common_asvs,
        ,
        drop = FALSE
    ]

    # ------------------------------------------------------------
    # Match samples to metadata
    # ------------------------------------------------------------

    sample_order <- intersect(
        metadata[["sample"]],
        colnames(asv)
    )

    if (length(sample_order) == 0) {
        stop("No samples were shared between ASV table and metadata.")
    }

    asv <- asv[
        ,
        sample_order,
        drop = FALSE
    ]

    metadata <- metadata[
        match(sample_order, metadata[["sample"]]),
        ,
        drop = FALSE
    ]

    # ------------------------------------------------------------
    # Calculate relative abundance
    # ------------------------------------------------------------

    relative_abundance <- sweep(
        asv,
        2,
        colSums(asv),
        "/"
    ) * 100

    # ------------------------------------------------------------
    # Convert abundance data to long format
    # ------------------------------------------------------------

    abundance_list <- list()

    for (i in seq_len(nrow(candidate_taxa))) {

        asv_id <- candidate_taxa[["ASV"]][i]

        if (!asv_id %in% rownames(relative_abundance)) {
            next
        }

        values <- as.numeric(
            relative_abundance[asv_id, ]
        )

        abundance_list[[length(abundance_list) + 1]] <- data.frame(
            ASV = asv_id,
            taxon = candidate_taxa[["taxon_label"]][i],
            sample = sample_order,
            treatment = metadata[["treatment"]],
            relative_abundance = values,
            stringsAsFactors = FALSE
        )
    }

    if (length(abundance_list) == 0) {
        stop("No candidate ASVs could be converted to abundance data.")
    }

    abundance_long <- do.call(
        rbind,
        abundance_list
    )

    # ------------------------------------------------------------
    # Write candidate tables
    # ------------------------------------------------------------

    write.table(
        candidate_taxa,
        file = "candidate_taxa.tsv",
        sep = "\\t",
        quote = FALSE,
        row.names = FALSE
    )

    write.table(
        abundance_long,
        file = "candidate_taxa_abundance.tsv",
        sep = "\\t",
        quote = FALSE,
        row.names = FALSE
    )

    # ------------------------------------------------------------
    # Prepare heatmap
    # ------------------------------------------------------------

    abundance_long[["taxon"]] <- factor(
        abundance_long[["taxon"]],
        levels = rev(candidate_taxa[["taxon_label"]])
    )

    abundance_long[["sample"]] <- factor(
        abundance_long[["sample"]],
        levels = sample_order
    )

    plot <- ggplot(
        abundance_long,
        aes(
            x = sample,
            y = taxon,
            fill = relative_abundance
        )
    ) +
        geom_tile() +
        facet_grid(
            ~ treatment,
            scales = "free_x",
            space = "free_x"
        ) +
        labs(
            title = "Exploratory candidate taxa",
            subtitle = "Relative abundance across soil-management treatments",
            x = "Sample",
            y = "Taxon",
            fill = "Relative abundance (%)"
        ) +
        theme_minimal(base_size = 11) +
        theme(
            axis.text.x = element_text(
                angle = 45,
                hjust = 1
            ),
            panel.grid = element_blank()
        )

    ggsave(
        "candidate_taxa_heatmap.png",
        plot = plot,
        width = 10,
        height = 6,
        dpi = 300
    )

    # ------------------------------------------------------------
    # Create summary
    # ------------------------------------------------------------

    beta_text <- readLines(
        beta_file,
        warn = FALSE
    )

    summary_lines <- c(
        "Exploratory candidate taxa analysis",
        "===================================",
        "",
        paste(
            "Exploratory candidates meeting |log2 fold change| >= 2:",
            nrow(candidates)
        ),
        paste(
            "Candidate taxa visualized:",
            nrow(candidate_taxa)
        ),
        "",
        "Candidates were ranked by absolute log2 fold change.",
        "The ten strongest exploratory candidates were retained for visualization.",
        "",
        "Candidate ASVs were linked to SILVA taxonomic assignments.",
        "Relative abundance was calculated for each candidate ASV across samples.",
        "",
        "Interpretation:",
        "These features are exploratory signals and should not be interpreted",
        "as statistically confirmed differential taxa.",
        "",
        "Community-level beta-diversity results:",
        beta_text
    )

    writeLines(
        summary_lines,
        con = "candidate_taxa_summary.txt"
    )

    RSCRIPT
    """
}