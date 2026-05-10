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

    candidates.sort(
        key=lambda member: ("/GCA_" not in member and "/GCF_" not in member, member)
    )
    return candidates[0]


def open_zip_files(zip_dir):
    paths = sorted(Path(zip_dir).glob("batch_*.zip"))
    if not paths:
        raise SystemExit(f"No batch zip files found in {zip_dir}")

    zip_files = []
    for path in paths:
        try:
            zip_files.append((path, zipfile.ZipFile(path)))
        except zipfile.BadZipFile:
            raise SystemExit(f"Invalid zip archive: {path}")

    return zip_files


accessions = read_accessions(snakemake.input.accessions)
output_dir = Path(snakemake.output.assemblies_dir)
manifest = Path(snakemake.output.manifest)

if output_dir.exists():
    shutil.rmtree(output_dir)
output_dir.mkdir(parents=True, exist_ok=True)
manifest.parent.mkdir(parents=True, exist_ok=True)

rows = []
missing = []
unpacked = 0
zip_files = open_zip_files(snakemake.input.zips)

try:
    for accession in accessions:
        output_fasta = output_dir / f"{accession}.fa"
        found_member = None
        found_zip = None
        found_zip_path = None

        for zip_path, zip_file in zip_files:
            member = find_fasta_member(zip_file, accession)
            if member:
                found_member = member
                found_zip = zip_file
                found_zip_path = zip_path
                break

        if found_member is None:
            missing.append(accession)
            rows.append(
                {
                    "assembly_ID": accession,
                    "fasta": "",
                    "zip_member": "",
                    "source_zip": "",
                    "status": "missing",
                }
            )
            continue

        with found_zip.open(found_member) as source, open(output_fasta, "wb") as target:
            shutil.copyfileobj(source, target)

        unpacked += 1
        rows.append(
            {
                "assembly_ID": accession,
                "fasta": str(output_fasta),
                "zip_member": found_member,
                "source_zip": str(found_zip_path),
                "status": "unpacked",
            }
        )
finally:
    for _, zip_file in zip_files:
        zip_file.close()

with open(manifest, "w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["assembly_ID", "fasta", "zip_member", "source_zip", "status"],
    )
    writer.writeheader()
    writer.writerows(rows)

if unpacked == 0:
    raise SystemExit("No assemblies were unpacked from NCBI batch zips.")

if missing:
    examples = ", ".join(missing[:5])
    raise SystemExit(f"Missing FASTA files for {len(missing)} accession(s): {examples}")
