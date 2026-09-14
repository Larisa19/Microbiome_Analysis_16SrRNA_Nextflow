nextflow.enable.dsl=2

process DADA2_FILTER {
tag "$sample_id"

publishDir 'results/dada2/filter', mode: 'copy'

input:

tuple val(sample_id), val(treatment), val(irrigation), val(cultivar), val(stage), path(read1), path(read2)

output:

tuple val(sample_id),
      val(treatment),
      val(irrigation),
      val(cultivar),
      val(stage),
      path("${sample_id}_F_filtered.fastq.gz"),
      path("${sample_id}_R_filtered.fastq.gz"),
      path("${sample_id}_filter_stats.tsv"),
      emit: filtered

script:

"""
Rscript - <<'EOF'

library(dada2)

fnFs <- "${read1}"
fnRs <- "${read2}"

filtFs <- "${sample_id}_F_filtered.fastq.gz"
filtRs <- "${sample_id}_R_filtered.fastq.gz"

out <- filterAndTrim(
    fnFs,
    filtFs,
    fnRs,
    filtRs,
    truncLen = c(220, 180),
    maxN = 0,
    maxEE = c(2, 2),
    truncQ = 2,
    rm.phix = TRUE,
    compress = TRUE,
    multithread = FALSE
)

write.table(
    out,
    file = "${sample_id}_filter_stats.tsv",
    sep = "\\t",
    quote = FALSE,
    col.names = TRUE
)

EOF
"""
}

process DADA2_LEARN_ERRORS {
tag "all_samples"

publishDir 'results/dada2/errors', mode: 'copy'

input:

path filtered_reads

output:

path "error_model_F.rds", emit: error_model_F
path "error_model_R.rds", emit: error_model_R

script:

"""
Rscript - <<'EOF'

library(dada2)

files <- list.files(
    ".",
    pattern = "_filtered.fastq.gz\$",
    full.names = TRUE
)

fnFs <- files[grepl("_F_filtered.fastq.gz\$", files)]
fnRs <- files[grepl("_R_filtered.fastq.gz\$", files)]

errF <- learnErrors(
    fnFs,
    multithread = TRUE,
    randomize = TRUE
)

errR <- learnErrors(
    fnRs,
    multithread = TRUE,
    randomize = TRUE
)

saveRDS(errF, "error_model_F.rds")
saveRDS(errR, "error_model_R.rds")

EOF
"""
    stub:
    """
    touch error_model_F.rds
    touch error_model_R.rds
    """
}
process DADA2_DENOISE {

    tag "$sample_id"

    publishDir 'results/dada2/denoise', mode: 'copy'

    input:

    tuple val(sample_id),
          val(treatment),
          val(irrigation),
          val(cultivar),
          val(stage),
          path(read1),
          path(read2)

    path error_model_F
    path error_model_R

    output:

    tuple val(sample_id),
          val(treatment),
          val(irrigation),
          val(cultivar),
          val(stage),
          path("${sample_id}_F_dada.rds"),
	  path("${sample_id}_R_dada.rds"),
	  path("${sample_id}_F_derep.rds"),
	  path("${sample_id}_R_derep.rds"),
emit: denoised

    script:

    """

    Rscript - <<'EOF'

    library(dada2)
	

    errF <- readRDS("${error_model_F}")
    errR <- readRDS("${error_model_R}")

    derepF <- derepFastq("${read1}")
    derepR <- derepFastq("${read2}")

    dadaF <- dada(
        "${read1}",
        err = errF,
        multithread = FALSE
    )

    dadaR <- dada(
        "${read2}",
        err = errR,
        multithread = FALSE
    )

    saveRDS(dadaF, "${sample_id}_F_dada.rds")
    saveRDS(dadaR, "${sample_id}_R_dada.rds")
    saveRDS(derepF, "${sample_id}_F_derep.rds")
    saveRDS(derepR, "${sample_id}_R_derep.rds")

EOF

    """
   stub:
   """
   touch ${sample_id}_F_dada.rds
   touch ${sample_id}_R_dada.rds
   touch ${sample_id}_F_derep.rds
   touch ${sample_id}_R_derep.rds
   """
}
process DADA2_MERGE {

    tag "$sample_id"

    publishDir 'results/dada2/merge', mode: 'copy'

    input:

    tuple val(sample_id),
          val(treatment),
          val(irrigation),
          val(cultivar),
          val(stage),
          path(dadaF),
	path(dadaR),
	path(derepF),
	path(derepR)

    output:

    tuple val(sample_id),
          val(treatment),
          val(irrigation),
          val(cultivar),
          val(stage),
          path("${sample_id}_merged.rds"),
          emit: merged

    script:

    """

    Rscript - <<'EOF'

    library(dada2)

   dadaF <- readRDS("${dadaF}")
   dadaR <- readRDS("${dadaR}")
   derepF <- readRDS("${derepF}")
   derepR <- readRDS("${derepR}")

mergers <- mergePairs(

    dadaF,
    derepF,
    dadaR,
    derepR
    )

    saveRDS(
        mergers,
        "${sample_id}_merged.rds"
    )

EOF

    """
    stub:
    """
    touch ${sample_id}_merged.rds
    """
}

process DADA2_TABLE {

    tag "ASV_table"

    publishDir 'results/dada2/table', mode: 'copy'

    input:

    path merged_files

    output:

    path "asv_table.rds", emit: asv_table
    path "asv_table.tsv", emit: asv_table_tsv

    script:

    """

    Rscript - <<'EOF'

    library(dada2)

    files <- list.files(
        ".",
        pattern = "_merged.rds\$",
        full.names = TRUE
    )

    sample_names <- sub(
        "_merged.rds\$",
        "",
        basename(files)
    )

    merged <- lapply(files, readRDS)
    names(merged) <- sample_names

    sequences <- unique(
        unlist(
            lapply(merged, function(x) x\$sequence)
        )
    )

    asv_table <- matrix(
        0,
        nrow = length(sequences),
        ncol = length(merged),
        dimnames = list(sequences, sample_names)
    )

    for (i in seq_along(merged)) {
        x <- merged[[i]]
        asv_table[x\$sequence, i] <- x\$abundance
    }

    saveRDS(
        asv_table,
        "asv_table.rds"
    )

    write.table(
        asv_table,
        file = "asv_table.tsv",
        sep = "\\t",
        quote = FALSE,
        col.names = NA
    )

    EOF

    """
    stub:
    """
    touch asv_table.rds
    touch asv_table.tsv
    """
}