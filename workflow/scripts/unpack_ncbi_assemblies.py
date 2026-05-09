#!/usr/bin/env python3
import argparse
import csv
import shutil
import zipfile
from pathlib import Path


FASTA_SUFFIXES = (".fa", ".fna", ".fasta")


def read_accessions(path):
    with open(path) as handle:
        return [line.strip() for line in handle if line.strip()]


def find_fasta_member(zip_file, accession):
    candidates = []
    accession_marker = f"/{accession}/"
    for member in zip_file.namelist():
        lower = member.lower()
        if not lower.endswith(FASTA_SUFFIXES):
            continue
        if accession_marker in member:
            candidates.append(member)
    if not candidates:
        return None
    candidates.sort(key=lambda member: ("/GCA_" not in member and "/GCF_" not in member, member))
    return candidates[0]


def main():
    parser = argparse.ArgumentParser(
        description="Unpack NCBI datasets genome zip into one FASTA per accession."
    )
    parser.add_argument("--zip", required=True, help="NCBI datasets zip")
    parser.add_argument("--accessions", required=True, help="Expected accession list")
    parser.add_argument("--output-dir", required=True, help="Output assembly FASTA directory")
    parser.add_argument("--manifest", required=True, help="Output unpack manifest CSV")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    manifest = Path(args.manifest)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest.parent.mkdir(parents=True, exist_ok=True)

    accessions = read_accessions(args.accessions)
    rows = []
    missing = []

    with zipfile.ZipFile(args.zip) as zip_file:
        for accession in accessions:
            member = find_fasta_member(zip_file, accession)
            output_fasta = output_dir / f"{accession}.fa"
            if member is None:
                missing.append(accession)
                rows.append(
                    {
                        "assembly_ID": accession,
                        "fasta": "",
                        "zip_member": "",
                        "status": "missing",
                    }
                )
                continue

            with zip_file.open(member) as source, open(output_fasta, "wb") as target:
                shutil.copyfileobj(source, target)

            rows.append(
                {
                    "assembly_ID": accession,
                    "fasta": str(output_fasta),
                    "zip_member": member,
                    "status": "unpacked",
                }
            )

    with open(manifest, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["assembly_ID", "fasta", "zip_member", "status"],
        )
        writer.writeheader()
        writer.writerows(rows)

    if missing:
        examples = ", ".join(missing[:5])
        raise SystemExit(f"Missing FASTA files for {len(missing)} accession(s): {examples}")


if __name__ == "__main__":
    main()
