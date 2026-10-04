# Vineyard Soil Microbial Ecology

## 16S rRNA amplicon analysis of vineyard soil microbial communities

### Project objective

This project investigates whether **soil management (tillage vs cover crop/no-tillage)** is associated with differences in microbial community composition in vineyard soil.

To isolate the effect of soil management, the analysis focuses on:

* **Cultivar:** Xynisteri
* **Stage:** harvest
* **Irrigation:** no irrigation
* **Treatments:** tillage vs cover crop/no-tillage
* **Samples:** 8 biological samples (4 per treatment)

The sequencing data are from **NCBI BioProject PRJEB40549**.

The project was implemented as a reproducible **Nextflow DSL2 workflow**, but the main focus is the biological analysis and the reasoning behind the analytical decisions.

---

# Biological question

> **Does soil management affect microbial community composition in non-irrigated Xynisteri vineyard soil at harvest?**

The analysis was deliberately restricted to samples sharing the same cultivar, developmental stage and irrigation condition. This reduces variation from other experimental factors and makes the comparison between the two soil-management treatments more focused.

---

# Samples

| Sample     | Treatment  | Irrigation    | Cultivar  | Stage   |
| ---------- | ---------- | ------------- | --------- | ------- |
| ERR4702541 | tillage    | no_irrigation | Xynisteri | harvest |
| ERR4702540 | tillage    | no_irrigation | Xynisteri | harvest |
| ERR4702539 | tillage    | no_irrigation | Xynisteri | harvest |
| ERR4702538 | tillage    | no_irrigation | Xynisteri | harvest |
| ERR4702533 | cover_crop | no_irrigation | Xynisteri | harvest |
| ERR4702532 | cover_crop | no_irrigation | Xynisteri | harvest |
| ERR4702531 | cover_crop | no_irrigation | Xynisteri | harvest |
| ERR4702530 | cover_crop | no_irrigation | Xynisteri | harvest |

Metadata:

```text
assets/samplesheet.csv
```

---

# Analysis workflow

The complete analysis follows this sequence:

```text
Raw paired-end FASTQ
        │
        ├── FastQC
        │
        └── Cutadapt
              │
              ▼
        Trimmed FASTQ
              │
              ├── FastQC
              │
              └── DADA2 filtering
                      │
                      ▼
                Filtered reads
                      │
             ┌────────┴─────────┐
             │                  │
             ▼                  ▼
       Error learning       Denoising
             │                  │
             └───────┬──────────┘
                     ▼
               Read merging
                     │
                     ▼
                 ASV table
                     │
        ┌────────────┼─────────────┐
        │            │             │
        ▼            ▼             ▼
      ASV QC    Alpha diversity  Taxonomy
        │            │             │
        │            │        Phylum abundance
        │            │             │
        │            │        Composition plot
        │
        └────── Rarefaction
                    │
                    ▼
             Rarefied ASV table
                    │
                    ▼
              Bray–Curtis
                    │
             ┌──────┼──────┐
             ▼      ▼      ▼
            PCoA  PERMANOVA PERMDISP

FastQC raw + FastQC trimmed
              │
              ▼
            MultiQC
```

---
## Step 1 — Raw read quality control

**Input:** Raw paired-end FASTQ files from `data/raw/`

**Tools:** FastQC + MultiQC

FastQC was used to assess the quality of the raw sequencing reads. MultiQC was then used to aggregate the FastQC results into an interactive report.

**Output:**

* FastQC reports: `results/qc/raw/fastqc/`
* MultiQC report: `results/qc/multiqc_raw/multiqc_raw_report.html`

### Result

Quality control was performed for all 8 samples before adapter trimming.

**[View interactive raw-read MultiQC report](https://larisa19.github.io/microbial_ecology/qc/raw/multiqc_raw_report.html)**

---

## Step 2 — Adapter trimming

**Input:** Raw paired-end FASTQ files from `data/raw/`

**Tool:** Cutadapt

Illumina/TruSeq adapter sequences were removed from the forward and reverse reads using Cutadapt.

**Output:** Trimmed paired-end FASTQ files in `results/cutadapt/`

---

## Step 3 — Trimmed read quality control

**Input:** Trimmed paired-end FASTQ files from `results/cutadapt/`

**Tools:** FastQC + MultiQC

FastQC was used to assess read quality after adapter removal. MultiQC was then used to aggregate the results into a separate interactive report.

**Output:**

* FastQC reports: `results/qc/trimmed/fastqc/`
* MultiQC report: `results/qc/multiqc_trimmed/multiqc_trimmed_report.html`

### Result

Quality control was performed for all 8 samples after adapter trimming.

**[View interactive trimmed-read MultiQC report](https://larisa19.github.io/microbial_ecology/qc/trimmed/multiqc_trimmed_report.html)**

## Step 4 — DADA2 quality filtering

**Input:** Trimmed paired-end FASTQ files

The trimmed reads were filtered to remove low-quality reads and reads containing ambiguous bases.

**Tool:** DADA2 `filterAndTrim()`

Parameters:

* `truncLen = c(220, 180)`
* `maxN = 0`
* `maxEE = c(2, 2)`
* `truncQ = 2`
* `rm.phix = TRUE`
* `compress = TRUE`

**Output:** Filtered forward and reverse FASTQ files and filtering statistics.

Files: `results/dada2/filter/`

---

## Step 5 — Learn the DADA2 error model

**Input:**

* Filtered forward reads from all samples
* Filtered reverse reads from all samples

The DADA2 error model was learned from the filtered reads to characterize sequencing errors.

**Tool:** DADA2 `learnErrors()`

Forward and reverse reads were modeled separately.

**Output:** `error_model_F.rds` and `error_model_R.rds`

Files: `results/dada2/errors/error_model_F.rds`, `results/dada2/errors/error_model_R.rds`

---

## Step 6 — Denoise reads

**Input:**

* Filtered forward and reverse reads
* Forward and reverse DADA2 error models

The filtered reads were denoised to distinguish biological sequence variation from sequencing errors.

**Tool:** DADA2 `dada()`

Forward and reverse reads were denoised separately using the corresponding error models.

**Output:** DADA2 denoising objects and dereplicated read objects for each sample.

Files: `results/dada2/denoise/`

---

## Step 7 — Paired-end merging

**Input:**

* Denoised forward reads
* Denoised reverse reads
* Corresponding dereplicated reads

The forward and reverse DADA2 results were merged to reconstruct the amplicon sequences.

**Tool:** DADA2 `mergePairs()`

**Output:** One merged object per sample.

Files: `results/dada2/merge/`

---

## Step 8 — Construct the ASV table

**Input:** Merged sequences from all samples

The merged sequences were combined into a sequence-by-sample abundance matrix.

**Tool:** DADA2

The workflow:

* reads all merged objects
* extracts unique ASV sequences
* creates a matrix with ASVs as rows and samples as columns
* fills the matrix with the corresponding sequence abundances

**Output:** `results/dada2/table/asv_table.rds`, `results/dada2/table/asv_table.tsv`

### Result

The final table contains:

* **8,953 ASVs**
* **625,399 total reads**

File: `results/dada2/table/asv_table.tsv`

The ASV table is the main input for the downstream ecological analyses.

---

## Step 9 — ASV quality control

**Input:**

* ASV table
* sample metadata

The ASV table was validated against the sample metadata before downstream ecological analyses.

**Tool:** Python

The workflow:

* reads the ASV table
* identifies all sample IDs
* counts total ASVs and reads per sample
* counts observed ASVs per sample
* reads the sample metadata
* compares sample IDs between the ASV table and metadata
* checks for invalid rows or missing samples

**Output:** 
`asv_qc_summary.tsv` — sample-level summary of sequencing reads and observed ASVs; `asv_qc_validation.txt` — validation report confirming that the ASV table and metadata are correctly aligned.

### Result

The validation confirmed:

* **8 samples** in the metadata
* **8 samples** in the ASV table
* **8 matching samples**
* **8,953 ASVs**
* **0 invalid rows**
* No samples missing from either the ASV table or metadata

Files: `results/asv_qc/asv_qc_summary.tsv`, `results/asv_qc/asv_qc_validation.txt`

The validated ASV table can now be used for downstream ecological analyses.

---
## Step 10 — Alpha diversity
Alpha diversity was calculated to describe microbial diversity within each sample.
**Input:**

* ASV table
* sample metadata

**Tool:** R

Two measures were calculated:

* **Observed ASVs** — number of detected ASVs per sample
* **Shannon diversity** — accounts for both richness and relative abundance

**Output:** `results/alpha_diversity/alpha_diversity.tsv`

### Result

Alpha-diversity values were calculated for all **8 samples**.

The results were used to compare within-sample diversity between the two soil-management treatments.

