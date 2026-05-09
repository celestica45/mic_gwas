#!/usr/bin/env python3
import argparse
import json
import subprocess
import time
from pathlib import Path

import pandas as pd
import requests


EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_ASSEMBLY_PAGE = "https://www.ncbi.nlm.nih.gov/assembly/{accession}/"


def ncbi_get(endpoint, params):
    response = requests.get(f"{EUTILS}/{endpoint}", params=params, timeout=60)
    response.raise_for_status()
    return response.json()


def find_assembly_id(accession):
    result = ncbi_get(
        "esearch.fcgi",
        {
            "db": "assembly",
            "term": accession,
            "retmode": "json",
        },
    )
    ids = result.get("esearchresult", {}).get("idlist", [])
    return ids[0] if ids else ""


def assembly_summary(accession):
    assembly_id = find_assembly_id(accession)
    if not assembly_id:
        return {}

    result = ncbi_get(
        "esummary.fcgi",
        {
            "db": "assembly",
            "id": assembly_id,
            "retmode": "json",
        },
    )
    return result.get("result", {}).get(assembly_id, {})


def get_value(summary, *keys):
    for key in keys:
        value = summary.get(key)
        if value:
            return value
    return ""


def get_stat(summary, key):
    stat = summary.get("assemblystats", {})
    return stat.get(key, "")


def datasets_summary(accession):
    try:
        result = subprocess.run(
            ["datasets", "summary", "genome", "accession", accession, "--as-json-lines"],
            check=True,
            capture_output=True,
            text=True,
            timeout=120,
        )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
        return {}

    for line in result.stdout.splitlines():
        line = line.strip()
        if line:
            return json.loads(line)
    return {}


def nested(data, *keys):
    current = data
    for key in keys:
        if not isinstance(current, dict):
            return ""
        current = current.get(key, "")
    return current


def biosample_attribute(report, name):
    attributes = nested(report, "assembly_info", "biosample", "attributes")
    if not isinstance(attributes, list):
        return ""
    for attribute in attributes:
        if attribute.get("name") == name:
            return attribute.get("value", "")
    return ""


def sample_id(report, label_or_db):
    sample_ids = nested(report, "assembly_info", "biosample", "sample_ids")
    if not isinstance(sample_ids, list):
        return ""
    for item in sample_ids:
        if item.get("label") == label_or_db or item.get("db") == label_or_db:
            return item.get("value", "")
    return ""


def bioproject_lineage(report):
    lineage = nested(report, "assembly_info", "bioproject_lineage")
    projects = []
    if not isinstance(lineage, list):
        return ""
    for group in lineage:
        for project in group.get("bioprojects", []):
            accession = project.get("accession", "")
            title = project.get("title", "")
            if accession or title:
                projects.append(f"{accession}:{title}")
    return ";".join(projects)


def main():
    parser = argparse.ArgumentParser(description="Fetch NCBI Assembly metadata.")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(args.metadata)
    rows = []

    for _, sample in metadata.iterrows():
        assembly_accession = sample.get("assembly_ID", "")
        report = datasets_summary(assembly_accession)
        summary = assembly_summary(assembly_accession)

        rows.append(
            {
                "isolate": sample.get("isolate", ""),
                "BioSample_ID": sample.get("BioSample_ID", ""),
                "assembly_ID": assembly_accession,
                "ncbi_assembly_accession": report.get("accession", "")
                or get_value(summary, "assemblyaccession"),
                "refseq_accession": report.get("paired_accession", ""),
                "assembly_name": nested(report, "assembly_info", "assembly_name")
                or get_value(summary, "assemblyname"),
                "organism": nested(report, "organism", "organism_name")
                or get_value(summary, "organism"),
                "tax_id": nested(report, "organism", "tax_id") or get_value(summary, "taxid"),
                "BioProject_ID": nested(report, "assembly_info", "bioproject_accession")
                or get_value(summary, "bioprojectaccn"),
                "BioProject_lineage": bioproject_lineage(report),
                "NCBI_BioSample": nested(report, "assembly_info", "biosample", "accession")
                or get_value(summary, "biosampleaccn"),
                "submitter": nested(report, "assembly_info", "submitter")
                or get_value(summary, "submitterorganization"),
                "strain": nested(report, "organism", "infraspecific_names", "strain")
                or biosample_attribute(report, "strain")
                or get_value(summary, "strain"),
                "sample_name": sample_id(report, "Sample name") or get_value(summary, "samplename"),
                "sra_accession": sample_id(report, "SRA") or get_value(summary, "sraaccession"),
                "biosample_description": nested(
                    report, "assembly_info", "biosample", "description", "title"
                ),
                "panel_id": biosample_attribute(report, "panel_id"),
                "assembly_level": nested(report, "assembly_info", "assembly_level")
                or get_value(summary, "assemblystatus"),
                "assembly_status": nested(report, "assembly_info", "assembly_status"),
                "assembly_type": nested(report, "assembly_info", "assembly_type"),
                "genome_representation": get_value(summary, "genomerepresentation"),
                "refseq_category": get_value(summary, "refseq_category"),
                "submission_date": nested(report, "assembly_info", "biosample", "submission_date")
                or get_value(summary, "submissiondate"),
                "release_date": nested(report, "assembly_info", "release_date")
                or get_value(summary, "releasedate"),
                "assembly_method": nested(report, "assembly_info", "assembly_method")
                or get_value(summary, "assemblymethod"),
                "genome_coverage": nested(report, "assembly_stats", "genome_coverage")
                or get_value(summary, "genomecoverage"),
                "sequencing_technology": nested(report, "assembly_info", "sequencing_tech")
                or get_value(summary, "sequencingtechnology"),
                "genome_size": nested(report, "assembly_stats", "total_sequence_length")
                or get_stat(summary, "total_sequence_length"),
                "total_ungapped_length": nested(
                    report, "assembly_stats", "total_ungapped_length"
                )
                or get_stat(summary, "total_ungapped_length"),
                "number_of_chromosomes": get_stat(summary, "number_of_chromosomes"),
                "number_of_scaffolds": nested(report, "assembly_stats", "number_of_scaffolds")
                or get_stat(summary, "number_of_scaffolds"),
                "scaffold_n50": nested(report, "assembly_stats", "scaffold_n50")
                or get_stat(summary, "scaffold_n50"),
                "number_of_contigs": nested(report, "assembly_stats", "number_of_contigs")
                or get_stat(summary, "number_of_contigs"),
                "contig_n50": nested(report, "assembly_stats", "contig_n50")
                or get_stat(summary, "contig_n50"),
                "gc_percent": nested(report, "assembly_stats", "gc_percent")
                or get_stat(summary, "gc_percent"),
                "ftp_path_genbank": get_value(summary, "ftppath_genbank"),
                "ftp_path_refseq": get_value(summary, "ftppath_refseq"),
            }
        )
        time.sleep(0.34)

    pd.DataFrame(rows).to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
