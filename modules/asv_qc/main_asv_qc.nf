process ASV_QC {

    tag "ASV table validation"

    input:
    path asv_table
    path samplesheet

    output:
    path "asv_qc_summary.tsv", emit: summary
    path "asv_qc_validation.txt", emit: validation

    script:

    """
    python3 - "$asv_table" "$samplesheet" <<'PY'

import sys
import csv

asv_file = sys.argv[1]
metadata_file = sys.argv[2]

# --------------------------------------------------
# Read ASV table
# --------------------------------------------------

with open(asv_file, newline="") as f:
    reader = csv.reader(f, delimiter="\t")

    header = next(reader)

    # First column contains ASV sequences
    sample_ids_table = header[1:]

    n_columns = len(header)

    total_asvs = 0
    sample_reads = [0] * len(sample_ids_table)
    sample_asvs = [0] * len(sample_ids_table)
    invalid_rows = 0

    for row in reader:

        if len(row) != n_columns:
            invalid_rows += 1
            continue

        total_asvs += 1

        for i, value in enumerate(row[1:]):

            try:
                count = int(value)
            except ValueError:
                count = 0

            sample_reads[i] += count

            if count > 0:
                sample_asvs[i] += 1


# --------------------------------------------------
# Read metadata
# --------------------------------------------------

with open(metadata_file, newline="") as f:
    reader = csv.DictReader(f)

    metadata_rows = list(reader)

sample_ids_metadata = [row["sample"] for row in metadata_rows]


# --------------------------------------------------
# Validate sample IDs
# --------------------------------------------------

table_set = set(sample_ids_table)
metadata_set = set(sample_ids_metadata)

missing_in_table = sorted(metadata_set - table_set)
missing_in_metadata = sorted(table_set - metadata_set)

matching_samples = sorted(table_set & metadata_set)


# --------------------------------------------------
# Write validation report
# --------------------------------------------------

with open("asv_qc_validation.txt", "w") as out:

    out.write("ASV TABLE VALIDATION\\n")
    out.write("====================\\n\\n")

    out.write(f"Samples in metadata: {len(metadata_set)}\\n")
    out.write(f"Samples in ASV table: {len(table_set)}\\n")
    out.write(f"Matching samples: {len(matching_samples)}\\n")
    out.write(f"Total ASVs: {total_asvs}\\n")
    out.write(f"Invalid rows: {invalid_rows}\\n\\n")

    out.write("Missing from ASV table:\\n")
    out.write(", ".join(missing_in_table) if missing_in_table else "None")
    out.write("\\n\\n")

    out.write("Missing from metadata:\\n")
    out.write(", ".join(missing_in_metadata) if missing_in_metadata else "None")
    out.write("\\n")


# --------------------------------------------------
# Write sample summary
# --------------------------------------------------

metadata_by_sample = {
    row["sample"]: row
    for row in metadata_rows
}

with open("asv_qc_summary.tsv", "w", newline="") as out:

    writer = csv.writer(out, delimiter="\t")

    writer.writerow([
        "sample",
        "treatment",
        "irrigation",
        "cultivar",
        "stage",
        "total_reads",
        "observed_ASVs"
    ])

    for i, sample in enumerate(sample_ids_table):

        metadata = metadata_by_sample.get(sample, {})

        writer.writerow([
            sample,
            metadata.get("treatment", "NA"),
            metadata.get("irrigation", "NA"),
            metadata.get("cultivar", "NA"),
            metadata.get("stage", "NA"),
            sample_reads[i],
            sample_asvs[i]
        ])

PY
    """

    stub:
    """
    touch asv_qc_summary.tsv
    touch asv_qc_validation.txt
    """
}