import shutil
import subprocess
import sys
import zipfile
from pathlib import Path


BATCH_SIZE = int(snakemake.params.get("batch_size", 10))
DATASETS_INCLUDE = snakemake.params.get("datasets_include", "genome")

output_dir = Path(snakemake.output.zips)
partial_dir = output_dir.with_name(output_dir.name + ".part")


def read_accessions(path):
    accessions = []
    for line in Path(path).read_text().splitlines():
        accession = line.strip()
        if accession:
            accessions.append(accession)
    return accessions


def assert_valid_zip(path):
    if not path.exists() or path.stat().st_size == 0:
        raise RuntimeError(f"NCBI Datasets did not create a usable archive: {path}")

    try:
        with zipfile.ZipFile(path) as handle:
            bad_member = handle.testzip()
    except zipfile.BadZipFile as error:
        raise RuntimeError(
            f"NCBI Datasets created an invalid zip archive: {path}"
        ) from error

    if bad_member is not None:
        raise RuntimeError(
            f"NCBI Datasets created a corrupt zip member in {path}: {bad_member}"
        )


legacy_zip = output_dir.with_name("ncbi_assemblies.zip")
legacy_partial_zip = output_dir.with_name("ncbi_assemblies.zip.part")
for path in [output_dir, partial_dir, legacy_zip, legacy_partial_zip]:
    if path.exists():
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()

partial_dir.mkdir(parents=True, exist_ok=True)
accessions = read_accessions(snakemake.input.accessions)

if not accessions:
    raise ValueError(f"No assembly accessions found in {snakemake.input.accessions}")

for batch_number, start in enumerate(range(0, len(accessions), BATCH_SIZE), start=1):
    batch = accessions[start : start + BATCH_SIZE]
    batch_file = partial_dir / f"batch_{batch_number:04d}_accessions.txt"
    batch_zip = partial_dir / f"batch_{batch_number:04d}.zip"
    batch_file.write_text("\n".join(batch) + "\n")

    command = [
        "datasets",
        "download",
        "genome",
        "accession",
        "--inputfile",
        str(batch_file),
        "--include",
        DATASETS_INCLUDE,
        "--filename",
        str(batch_zip),
    ]

    print(
        f"Downloading assembly batch {batch_number} "
        f"({start + 1}-{start + len(batch)} of {len(accessions)})",
        file=sys.stderr,
    )

    completed = subprocess.run(command, text=True)
    if completed.returncode != 0:
        raise subprocess.CalledProcessError(completed.returncode, command)

    assert_valid_zip(batch_zip)

partial_dir.replace(output_dir)
