# Soil Microbial Ecology: Tillage vs Cover Crop

> **Work in progress** — This project is currently under development.

## Overview

This project develops a reproducible **Nextflow DSL2 workflow** for the analysis of soil microbial communities from a vineyard system using **16S rRNA amplicon sequencing data**.

The main biological question is:

> **Does soil management (tillage vs cover crop/no-tillage) affect the microbial community composition of non-irrigated Xynisteri vineyard soil at harvest?**

The workflow connects reproducible bioinformatics processing with downstream ecological and taxonomic analysis.

## Study Design

The current analysis includes 8 soil samples:

* 4 samples: tillage
* 4 samples: cover crop/no-tillage

All selected samples share:

* Cultivar: Xynisteri
* Sampling stage: harvest
* Irrigation: no irrigation

This design reduces variation from cultivar, sampling stage, and irrigation when comparing soil management treatments.

## Dataset

Sequencing data are paired-end Illumina MiSeq reads obtained from the public **NCBI/ENA repositories**.

Raw FASTQ files are intentionally **not included** in this repository.

Sample metadata and sequencing accessions are provided in:

```text
assets/samplesheet.csv
```

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
              |                                  v
              |                         Rarefaction
              |                                  |
              +----------------------------------+
                                             |
                                             v
                                      Beta diversity
                                             |
                                             v
                                  Bray-Curtis + PCoA
                                             |
                                             v
                              PERMANOVA + PERMDISP
                                             |
                                             v
                                  Taxonomic assignment
                                             |
                                             v
                                  Phylum composition
                                             |
                                             v
                              Treatment-associated taxa
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
* Rarefied ASV table
* Bray-Curtis distance matrix
* PCoA ordination
* PERMANOVA
* PERMDISP
* SILVA-based taxonomic assignment
* Phylum-level relative abundance
* Phylum composition visualization
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

Alpha diversity was calculated for each sample using:

* **Observed ASVs** — a measure of microbial richness
* **Shannon diversity** — a measure incorporating both richness and relative abundance distribution

Because sequencing depth varies between samples, a standardized sequencing depth was used for downstream diversity analysis through rarefaction.

The rarefied ASV table is used for the beta-diversity analysis.

## Beta Diversity and Community Structure

Community-level differences were assessed using **Bray-Curtis dissimilarity** based on the rarefied ASV table.

Ordination was performed using **Principal Coordinates Analysis (PCoA)**.

The first two PCoA axes explained:

* PCoA1: 18.71%
* PCoA2: 18.41%

### PERMANOVA

A PERMANOVA test with 999 permutations was used to evaluate whether microbial community composition differed between soil-management treatments.

The current analysis produced:

* R² = 0.170
* F = 1.233
* p = 0.03

This indicates that treatment explained approximately **17% of the variation in Bray-Curtis community composition** in this dataset.

### PERMDISP

PERMDISP was used to assess whether differences in within-group dispersion could account for the PERMANOVA result.

Current result:

* F = 0.415
* p = 0.552

There is no evidence in this dataset of a significant difference in within-group dispersion between treatments.

Because the current comparison contains only four samples per treatment, these results are considered **exploratory** and should be interpreted cautiously.

## Taxonomic Assignment

ASVs were taxonomically assigned using the **SILVA 138.1 reference database** with the DADA2 `assignTaxonomy` method.

The resulting taxonomy table contains assignments at:

* Kingdom
* Phylum
* Class
* Order
* Family
* Genus

The SILVA reference database is not stored in the GitHub repository because of its file size. It is downloaded separately and excluded through `.gitignore`.

## Phylum-Level Community Composition

Relative abundance was calculated at the phylum level by aggregating ASV abundances according to their taxonomic assignment.

The current dataset is dominated by several major bacterial and archaeal groups, including:

* Actinobacteriota
* Proteobacteria
* Acidobacteriota
* Chloroflexi
* Bacteroidota
* Crenarchaeota
* Verrucomicrobiota
* Myxococcota
* Firmicutes

A treatment-level composition plot summarizes the relative abundance of the most abundant phyla across tillage and cover crop/no-tillage treatments.

## Planned Analysis

The next stage will focus on identifying **treatment-associated taxa** rather than adding additional general preprocessing steps.

Planned analyses include:

### 1. Taxonomic treatment comparison

Compare taxonomic composition between:

* Tillage
* Cover crop/no-tillage

at appropriate taxonomic levels.

### 2. Treatment-associated taxa

Identify microbial groups showing differences in relative abundance between soil-management treatments.

### 3. Ecological interpretation

Interpret the observed community-level and taxonomic patterns in the context of:

* Soil management
* Vineyard soil ecology
* Microbial diversity
* Potential effects of cover cropping and reduced soil disturbance

The final interpretation will distinguish exploratory statistical associations from biological conclusions.

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
* SILVA 138.1 training reference for taxonomic assignment

## Reproducibility

The workflow is implemented using **Nextflow DSL2** and organized into modular processes.

Raw sequencing data and large reference databases are not stored in the repository.

The expected input structure is:

```text
data/raw/
assets/samplesheet.csv
```

The complete workflow can be launched with:

```bash
nextflow run main.nf
```

For a previously cached workflow run:

```bash
nextflow run main.nf -resume
```

Individual downstream modules can also be executed independently using existing results. This allows analysis and visualization steps to be reproduced without re-running the entire upstream pipeline.

## Project Structure

```text
microbial_ecology/
├── assets/
│   └── samplesheet.csv
├── data/
│   └── raw/                              # Raw FASTQ files (not tracked)
├── modules/
│   ├── dada2/
│   │   └── main_dada2.nf
│   ├── asv_qc/
│   │   └── main_asv_qc.nf
│   ├── alpha_diversity/
│   │   └── main_alpha_diversity.nf
│   ├── rarefaction/
│   │   └── main_rarefaction.nf
│   ├── beta_diversity/
│   │   └── main_beta_diversity.nf
│   ├── figures/
│   │   └── main_figures.nf
│   └── taxonomy/
│       ├── main_taxonomy.nf
│       ├── main_taxonomic_abundance.nf
│       └── main_taxonomy_plot.nf
├── docs/
│   ├── analysis_notes.md
│   └── table_asv/
│       └── asv_table.tsv
├── results/
│   ├── asv_qc/
│   ├── alpha_diversity/
│   ├── beta_diversity/
│   ├── rarefaction/
│   ├── figures/
│   └── taxonomy/
├── main.nf
├── nextflow.config
├── README.md
├── PROJECT_SUMMARY.md
└── .gitignore
```

## Next Steps

The core processing and community-level analysis are now in place.

The next analytical stage is:

**Taxonomic comparison → Treatment-associated taxa → Ecological interpretation**

The final goal is to build a reproducible analysis connecting **soil management practices with vineyard soil microbial community structure, diversity, and taxonomic composition**.
