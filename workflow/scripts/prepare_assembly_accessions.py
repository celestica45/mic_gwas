#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Write unique assembly accessions from cleaned metadata."
    )
    parser.add_argument("--metadata", required=True, help="Clean metadata with assembly_ID")
    parser.add_argument("--output", required=True, help="Assembly accession list")
    parser.add_argument("--assembly-column", default="assembly_ID")
    args = parser.parse_args()

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(args.metadata)
    if args.assembly_column not in metadata.columns:
        raise SystemExit(f"Missing required column: {args.assembly_column}")

    accessions = (
        metadata[args.assembly_column]
        .dropna()
        .astype(str)
        .str.strip()
    )
    accessions = sorted(accession for accession in accessions.unique() if accession)
    if not accessions:
        raise SystemExit("No assembly accessions found.")

    with open(args.output, "w") as handle:
        for accession in accessions:
            handle.write(f"{accession}\n")


if __name__ == "__main__":
    main()
