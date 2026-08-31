# Soil Microbial Ecology: Tillage vs Cover Crop

> **Work in progress** — This project is currently under development.

## Overview

This project develops a reproducible **Nextflow workflow** for the analysis of
soil microbial communities from a vineyard system using **16S rRNA amplicon
sequencing data**.

The main biological question is:

> Does soil management (tillage vs cover crop/no-tillage) affect the
> microbial community composition of non-irrigated Xynisteri vineyard soil
> at harvest?

## Study Design

The current analysis includes 8 soil samples:

* 4 samples: tillage
* 4 samples: cover crop/no-tillage

All selected samples share:

* Cultivar: Xynisteri
* Sampling stage: harvest
* Irrigation: no irrigation

## Dataset

Sequencing data are paired-end Illumina MiSeq reads obtained from the
NCBI/ENA public repositories.

Raw FASTQ files are intentionally **not included** in this repository.

Sample metadata and sequencing accessions are provided in the samplesheet.

## Workflow

The current Nextflow workflow includes:

```text
   Raw FASTQ
   |
   +----------------------+
   |                      |
   v                      |
FastQC                    |
   |                      |
   v                      |
Cutadapt                  |
   |                      |
   v                      |
FastQC after trimming     |
   |                      |
   +----------+-----------+
              |
              v
           MultiQC


Raw FASTQ
   |
   v
DADA2 filtering
   |
   v
Error model learning
   |
   v
DADA2 denoising
   |
   v
Paired-end merging
   |
   v
ASV table

```

## Current Status

The preprocessing and initial DADA2 workflow have been implemented and
successfully tested on the current 8-sample dataset.

The current dataset contains:

* 8 samples
* 8,953 inferred ASVs
* 625,399 total merged reads

Current outputs include:

* Raw-read quality control
* Adapter trimming
* Post-trimming quality control
* DADA2 filtering
* DADA2 error-model estimation
* Denoising
* Paired-end read merging
* ASV abundance table
* MultiQC quality control report: [`multiqc_report.html`](results/qc/multiqc/multiqc_report.html)

Downstream microbial community analysis is currently in progress.

Planned analyses include:

* ASV quality assessment
* Taxonomic assignment
* Alpha diversity
* Beta diversity
* Community composition
* Comparison of tillage vs cover crop treatments

## Requirements

The workflow requires:

- Nextflow
- Java
- FastQC
- Cutadapt
- MultiQC
- R with the DADA2 package

The workflow is designed to run locally and can be adapted to HPC environments.

## Reproducibility

The workflow is implemented using **Nextflow DSL2** and is designed to
separate workflow orchestration from individual analysis modules.

Raw sequencing data are not included in the repository. The workflow
expects paired-end FASTQ files in `data/raw/` and sample metadata in
`assets/samplesheet.csv`.

To run the workflow:

```bash
nextflow run main.nf

## Project Structure

```text
microbial_ecology/
├── assets/
│   └── samplesheet.csv
├── data/
│   └── raw/                  # Raw FASTQ files (not tracked)
├── modules/
│   └── dada2/
│       └── main_dada2.nf
├── results/
│   └── dada2/
├── main.nf
├── README.md
├── PROJECT_SUMMARY.md
└── .gitignore
```

## Next Steps

The next stage of the project will focus on downstream analysis of the
ASV table, including taxonomic classification and community-level
comparisons between soil management treatments.

