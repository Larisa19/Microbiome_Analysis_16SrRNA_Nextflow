# Soil Microbial Ecology: Tillage vs Cover Crop

> **Work in progress** — This project is currently under development.

## Overview

This project develops a reproducible **Nextflow workflow** for the analysis of soil microbial communities from a vineyard system using **16S rRNA amplicon sequencing data**.

The main biological question is:

> **Does soil management (tillage vs cover crop/no-tillage) affect the microbial community composition of non-irrigated Xynisteri vineyard soil at harvest?**

The workflow is designed to connect reproducible bioinformatics processing with downstream ecological analysis.

## Study Design

The current analysis includes 8 soil samples:

* 4 samples: tillage
* 4 samples: cover crop/no-tillage

All selected samples share:

* Cultivar: Xynisteri
* Sampling stage: harvest
* Irrigation: no irrigation

This design was selected to reduce variation from cultivar, sampling stage, and irrigation when comparing soil management treatments.

## Dataset

Sequencing data are paired-end Illumina MiSeq reads obtained from the public **NCBI/ENA repositories**.

Raw FASTQ files are intentionally **not included** in this repository.

Sample metadata and sequencing accessions are provided in:

`assets/samplesheet.csv`

## Workflow

```text
                    Raw paired-end FASTQ
                            |
              +-------------+-------------+
              |                           |
              v                           v
           FastQC                       DADA2
              |                           |
              v                           v
          Cutadapt                 Read filtering
              |                           |
              v                           v
     FastQC after trimming          Error modeling
              |                           |
              |                           v
              |                       Denoising
              |                           |
              |                           v
              |                  Paired-end merging
              |                           |
              |                           v
              |                       ASV table
              |                           |
              |                    +------+------+
              |                    |             |
              |                    v             v
              |                 ASV QC     Alpha diversity
              |                                  |
              +------------------+---------------+
                                 |
                                 v
                               MultiQC
                                 |
                                 v
                    Standardized sequencing depth
                                 |
                                 v
                         Beta diversity
                                 |
                                 v
                     Community composition
                                 |
                                 v
                       Taxonomic assignment
                                 |
                                 v
                       Treatment comparison
                                 |
                                 v
                    Ecological interpretation
```

The workflow is implemented using **Nextflow DSL2**, with individual analysis steps organized into reusable modules.

## Current Status

The preprocessing and initial ecological analysis have been successfully implemented and tested on the current 8-sample dataset.

### Sequencing and ASV processing

* 8 samples
* 8,953 inferred ASVs
* 625,399 total reads
* Paired-end Illumina MiSeq data
* DADA2-based ASV inference

The workflow currently produces:

* Raw-read quality control
* Adapter trimming
* Post-trimming quality control
* DADA2 read filtering
* Error-model estimation
* Denoising
* Paired-end merging
* ASV abundance table
* ASV table validation
* Observed ASV richness
* Shannon diversity
* MultiQC report

### ASV validation

The ASV table is automatically checked against the sample metadata.

Current validation:

* Samples in metadata: 8
* Samples in ASV table: 8
* Matching samples: 8
* Invalid rows: 0
* Missing samples: none

## Alpha Diversity

Alpha diversity is calculated for each sample using:

* **Observed ASVs** — a measure of microbial richness
* **Shannon diversity** — a measure incorporating both richness and relative abundance distribution

Sequencing depth varies between samples, so alpha-diversity comparisons will be performed after **standardizing sequencing depth**.

This step is important for distinguishing biological differences from differences caused by sequencing effort.

## Planned Ecological Analysis

The downstream analysis follows the biological question from diversity to community composition:

### 1. Standardized alpha diversity

Compare microbial richness and diversity between:

* Tillage
* Cover crop/no-tillage

after controlling for sequencing depth.

### 2. Beta diversity

Assess differences in overall microbial community structure between soil-management treatments using:

* Community distance matrices
* Ordination
* Statistical comparison of treatment groups

### 3. Community composition

Determine whether soil-management treatment is associated with systematic differences in microbial community composition.

### 4. Taxonomic assignment

Assign taxonomy to ASVs and characterize the bacterial groups present in the vineyard soil communities.

### 5. Treatment-associated taxa

Identify taxa or microbial groups associated with differences between tillage and cover crop/no-tillage treatments.

### 6. Ecological interpretation

Connect the observed microbial patterns to the original ecological question:

> **Does soil management affect the structure and diversity of the vineyard soil microbial community?**

The analysis will distinguish sequencing-depth effects from biological differences and interpret the results in the context of soil microbial ecology and vineyard management.

## Quality Control

Quality-control results are summarized using **MultiQC**.

The current MultiQC report is available through GitHub Pages:

`https://larisa19.github.io/microbial_ecology/multiqc/`

## Requirements

The workflow requires:

* Nextflow
* Java
* FastQC
* Cutadapt
* MultiQC
* R
* DADA2

## Reproducibility

The workflow is implemented using **Nextflow DSL2** and organized into modular processes.

Raw sequencing data are not stored in the repository.

The expected input structure is:

```text
data/raw/
assets/samplesheet.csv
```

Run the complete workflow with:

```bash
nextflow run main.nf
```

To resume a previous run using cached results:

```bash
nextflow run main.nf -resume
```

## Project Structure

```text
microbial_ecology/
├── assets/
│   └── samplesheet.csv
├── data/
│   └── raw/                         # Raw FASTQ files (not tracked)
├── modules/
│   ├── dada2/
│   │   └── main_dada2.nf
│   ├── asv_qc/
│   │   └── main_asv_qc.nf
│   └── alpha_diversity/
│       └── main_alpha_diversity.nf
├── docs/
│   ├── analysis_notes.md
│   └── table_asv/
│       └── asv_table.tsv
├── results/
│   ├── asv_qc_summary.tsv
│   ├── asv_qc_validation.txt
│   └── alpha_diversity.tsv
├── main.nf
├── nextflow.config
├── README.md
├── PROJECT_SUMMARY.md
└── .gitignore
```

## Next Steps

The next stage of the workflow will standardize sequencing depth for alpha-diversity analysis.

This will be followed by:

**Alpha diversity → Beta diversity → Community composition → Taxonomic assignment → Treatment-associated taxa → Ecological interpretation**

The final goal is to build a reproducible analysis connecting **soil management practices with vineyard soil microbial community structure and diversity**.
