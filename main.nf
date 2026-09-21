nextflow.enable.dsl=2

include { DADA2_FILTER; DADA2_LEARN_ERRORS; DADA2_DENOISE; DADA2_MERGE; DADA2_TABLE } from './modules/dada2/main_dada2.nf'
include { ASV_QC } from './modules/asv_qc/main_asv_qc.nf'
include { ALPHA_DIVERSITY } from './modules/alpha_diversity/main_alpha_diversity.nf'

process FASTQC_RAW {
tag "$sample_id"

input:

tuple val(sample_id), val(treatment), val(irrigation), val(cultivar), val(stage), path(reads)

output:

path "*_fastqc.*", emit: fastqc_raw

script:

"""
fastqc ${reads}
"""
}

process CUTADAPT {

tag "$sample_id"

input:

tuple val(sample_id), val(treatment), val(irrigation), val(cultivar), val(stage), path(reads)

output:

tuple val(sample_id),
      val(treatment),
      val(irrigation),
      val(cultivar),
      val(stage),
      path("${sample_id}_trimmed_1.fastq"),
      path("${sample_id}_trimmed_2.fastq"),
      emit: trimmed_reads

script:

"""
cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -o ${sample_id}_trimmed_1.fastq \
    -p ${sample_id}_trimmed_2.fastq \
    ${reads[0]} ${reads[1]}
"""
}

process FASTQC_TRIMMED {
tag "$sample_id"

input:

tuple val(sample_id),
      val(treatment),
      val(irrigation),
      val(cultivar),
      val(stage),
      path(read1),
      path(read2)

output:

path "*_fastqc.*", emit: fastqc_trimmed

script:

"""
fastqc ${read1} ${read2}
"""
}

process MULTIQC {
input:

path fastqc_results

output:

path "multiqc_report.html"

script:

"""
multiqc . -o . -n multiqc_report.html
"""
}

workflow {
samples = Channel
    .fromPath('assets/samplesheet.csv')
    .splitCsv(header: true)
    .map { row ->

        def r1 = file("data/raw/${row.err}_1.fastq")
        def r2 = file("data/raw/${row.err}_2.fastq")

        tuple(
            row.sample,
            row.treatment,
            row.irrigation,
            row.cultivar,
            row.stage,
            [r1, r2]
        )
    }

FASTQC_RAW(samples)

CUTADAPT(samples)

FASTQC_TRIMMED(CUTADAPT.out.trimmed_reads)

DADA2_FILTER(CUTADAPT.out.trimmed_reads)

filtered_for_denoise = DADA2_FILTER.out.filtered
    .map { sample_id, treatment, irrigation, cultivar, stage, read1, read2, stats ->
        tuple(sample_id, treatment, irrigation, cultivar, stage, read1, read2)
    }

filtered_for_errors = DADA2_FILTER.out.filtered
    .map { sample_id, treatment, irrigation, cultivar, stage, read1, read2, stats ->
        tuple(read1, read2)
    }
    .collect()


DADA2_LEARN_ERRORS(filtered_for_errors)

error_model_F = DADA2_LEARN_ERRORS.out.error_model_F.first()
error_model_R = DADA2_LEARN_ERRORS.out.error_model_R.first()

DADA2_DENOISE(
    filtered_for_denoise,
    error_model_F,
    error_model_R
)

DADA2_MERGE(DADA2_DENOISE.out.denoised)

merged_files = DADA2_MERGE.out.merged
    .map { sample_id, treatment, irrigation, cultivar, stage, merged ->
        merged
    }
    .collect()

DADA2_TABLE(merged_files)

ASV_QC(
    DADA2_TABLE.out.asv_table_tsv,
    Channel.value(file('assets/samplesheet.csv'))
)
ALPHA_DIVERSITY(
    DADA2_TABLE.out.asv_table_tsv,
    Channel.value(file('assets/samplesheet.csv'))
)
all_fastqc = FASTQC_RAW.out.fastqc_raw
    .mix(FASTQC_TRIMMED.out.fastqc_trimmed)
    .collect()

MULTIQC(all_fastqc)

}


