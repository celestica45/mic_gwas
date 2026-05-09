#!/usr/bin/env python3
import argparse
import csv


def read_phenotype_samples(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "wgs_id" not in (reader.fieldnames or []):
            raise SystemExit("Phenotype table is missing required column: wgs_id")
        return {row["wgs_id"].strip() for row in reader if row.get("wgs_id", "").strip()}


def read_rtab_samples(path):
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise SystemExit(f"Empty Rtab file: {path}")
    if not header or header[0] != "Gene":
        raise SystemExit("Rtab header must start with Gene")
    return {value.strip() for value in header[1:] if value.strip()}


def read_distance_samples(path):
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise SystemExit(f"Empty distance matrix: {path}")
    return {value.strip() for value in header[1:] if value.strip()}


def read_metadata_header(path):
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            return next(reader)
        except StopIteration:
            raise SystemExit(f"Empty metadata table: {path}")


def main():
    parser = argparse.ArgumentParser(description="Validate pyseer gene GWAS inputs.")
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--rtab", required=True)
    parser.add_argument("--distances", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--covariates", nargs="*", default=[])
    args = parser.parse_args()

    phenotype_samples = read_phenotype_samples(args.phenotype)
    rtab_samples = read_rtab_samples(args.rtab)
    distance_samples = read_distance_samples(args.distances)

    shared = phenotype_samples & rtab_samples & distance_samples
    if not shared:
        raise SystemExit(
            "No shared samples across phenotype, Rtab, and Mash distance matrix. "
            "Check sample labels."
        )

    missing_distance = sorted((phenotype_samples & rtab_samples) - distance_samples)
    if missing_distance:
        examples = ", ".join(missing_distance[:10])
        raise SystemExit(
            "Samples present in phenotype/Rtab but missing from Mash distance matrix: "
            f"{examples}"
        )

    metadata_header = read_metadata_header(args.metadata)
    missing_covariates = [
        covariate for covariate in args.covariates if covariate not in metadata_header
    ]
    if missing_covariates:
        raise SystemExit(
            "Missing configured covariate column(s): "
            + ", ".join(missing_covariates)
        )


if __name__ == "__main__":
    main()
