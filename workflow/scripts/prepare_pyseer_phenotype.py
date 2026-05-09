#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def read_metadata(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        columns = set(reader.fieldnames or [])
    return rows, columns


def main():
    parser = argparse.ArgumentParser(
        description="Create GWAS phenotype, covariate, and sample files from cleaned metadata."
    )
    parser.add_argument("--metadata", required=True, help="Clean metadata with assemblies CSV")
    parser.add_argument("--phenotype-tsv", required=True, help="Output GitHub-style phenotype TSV")
    parser.add_argument("--metadata-tsv", required=True, help="Output GitHub-style metadata/covariate TSV")
    parser.add_argument("--sample-column", default="assembly_ID")
    parser.add_argument("--log2-column", default="mic_log2")
    parser.add_argument("--log2-column-name", required=True)
    parser.add_argument("--covariates", nargs="*", default=[])
    args = parser.parse_args()

    for output in [args.phenotype_tsv, args.metadata_tsv]:
        Path(output).parent.mkdir(parents=True, exist_ok=True)

    rows, columns = read_metadata(args.metadata)
    required = {args.sample_column, args.log2_column}
    missing = required - columns
    if missing:
        raise SystemExit(f"Missing required column(s): {', '.join(sorted(missing))}")

    missing_covariates = sorted(set(args.covariates) - columns)
    if missing_covariates:
        raise SystemExit(
            f"Missing configured covariate column(s): {', '.join(missing_covariates)}"
        )

    records = {}
    conflicts = []
    for row in rows:
        sample = row.get(args.sample_column, "").strip()
        log2 = row.get(args.log2_column, "").strip()

        if not sample or not log2:
            continue

        try:
            float(log2)
        except ValueError:
            continue

        if sample in records:
            previous = records[sample]
            previous_log2 = previous.get(args.log2_column, "").strip()
            if log2 != previous_log2:
                conflicts.append(sample)
            continue

        records[sample] = row

    if conflicts:
        examples = ", ".join(sorted(set(conflicts))[:5])
        raise SystemExit(
            "Some sample labels have conflicting log2 phenotype values. "
            f"Examples: {examples}"
        )

    ordered_samples = sorted(records)

    with open(args.phenotype_tsv, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["wgs_id", args.log2_column_name])
        for sample in ordered_samples:
            row = records[sample]
            writer.writerow([sample, row.get(args.log2_column, "")])

    with open(args.metadata_tsv, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["wgs_id"] + args.covariates)
        for sample in ordered_samples:
            row = records[sample]
            writer.writerow([sample] + [row.get(column, "") for column in args.covariates])


if __name__ == "__main__":
    main()
