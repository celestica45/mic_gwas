#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def base_accession(value):
    return value.rsplit(".", 1)[0]


def build_label_map(manifest_path):
    label_map = {}
    duplicates = []

    with open(manifest_path, newline="") as handle:
        reader = csv.DictReader(handle)
        if "assembly_ID" not in (reader.fieldnames or []):
            raise SystemExit("Manifest is missing required column: assembly_ID")

        for row in reader:
            assembly_id = row["assembly_ID"].strip()
            if not assembly_id:
                continue
            base_id = base_accession(assembly_id)
            if base_id in label_map and label_map[base_id] != assembly_id:
                duplicates.append(base_id)
            label_map[base_id] = assembly_id

    if duplicates:
        examples = ", ".join(sorted(set(duplicates))[:5])
        raise SystemExit(f"Duplicate base assembly IDs in manifest: {examples}")

    return label_map


def normalize_label(label, label_map):
    if label in label_map.values():
        return label
    if label in label_map:
        return label_map[label]
    raise KeyError(label)


def main():
    parser = argparse.ArgumentParser(
        description="Normalize Mash distance matrix labels to versioned assembly IDs."
    )
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    label_map = build_label_map(args.manifest)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    missing = set()
    with open(args.input, newline="") as in_handle, open(
        output, "w", newline=""
    ) as out_handle:
        reader = csv.reader(in_handle, delimiter="\t")
        writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")

        try:
            header = next(reader)
        except StopIteration:
            raise SystemExit(f"Empty Mash distance matrix: {args.input}")

        normalized_header = [header[0]]
        for label in header[1:]:
            try:
                normalized_header.append(normalize_label(label, label_map))
            except KeyError:
                missing.add(label)
                normalized_header.append(label)

        writer.writerow(normalized_header)

        for row in reader:
            if not row:
                continue
            try:
                row[0] = normalize_label(row[0], label_map)
            except KeyError:
                missing.add(row[0])
            writer.writerow(row)

    if missing:
        examples = ", ".join(sorted(missing)[:10])
        output.unlink(missing_ok=True)
        raise SystemExit(f"Mash labels not found in assembly manifest: {examples}")


if __name__ == "__main__":
    main()
