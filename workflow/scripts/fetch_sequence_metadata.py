#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd
import requests


ENA_FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"


ENA_FIELDS = [
    "study_accession",
    "secondary_study_accession",
    "sample_accession",
    "secondary_sample_accession",
    "experiment_accession",
    "run_accession",
    "submission_accession",
    "scientific_name",
    "tax_id",
    "library_layout",
    "library_name",
    "library_strategy",
    "library_source",
    "library_selection",
    "instrument_platform",
    "instrument_model",
    "read_count",
    "base_count",
    "fastq_ftp",
    "fastq_md5",
    "submitted_ftp",
    "sample_alias",
    "sample_description",
    "center_name",
    "country",
    "collection_date",
    "host",
    "isolation_source",
]


def parse_tsv(text):
    if not text.strip():
        return pd.DataFrame()

    rows = [line.split("\t") for line in text.strip().splitlines()]
    if len(rows) <= 1:
        return pd.DataFrame()

    return pd.DataFrame(rows[1:], columns=rows[0])


def ena_report(accession):
    if not accession or pd.isna(accession):
        return pd.DataFrame()

    response = requests.get(
        ENA_FILEREPORT,
        params={
            "accession": accession,
            "result": "read_run",
            "fields": ",".join(ENA_FIELDS),
            "format": "tsv",
        },
        timeout=60,
    )
    if not response.ok:
        return pd.DataFrame()

    return parse_tsv(response.text)


def query_by_isolate(isolate):
    if not isolate or pd.isna(isolate):
        return pd.DataFrame()

    response = requests.get(
        "https://www.ebi.ac.uk/ena/portal/api/search",
        params={
            "result": "read_run",
            "query": f'sample_alias="{isolate}"',
            "fields": ",".join(ENA_FIELDS),
            "format": "tsv",
        },
        timeout=60,
    )
    if not response.ok:
        return pd.DataFrame()

    return parse_tsv(response.text)


def split_fastqs(value):
    if not value or pd.isna(value):
        return []
    return [item.strip() for item in str(value).split(";") if item.strip()]


def to_https(link):
    if not link:
        return ""
    if link.startswith("ftp://") or link.startswith("http"):
        return link.replace("ftp://", "https://")
    return "https://" + link


def ncbi_link(accession):
    if not accession or pd.isna(accession):
        return ""
    accession = str(accession)
    if accession.startswith("PRJ"):
        return f"https://www.ncbi.nlm.nih.gov/bioproject/{accession}"
    if accession.startswith("SRP"):
        return f"https://www.ncbi.nlm.nih.gov/sra?term={accession}"
    return ""


def pubmed_links(pubmed_ids):
    links = []
    for pubmed_id in str(pubmed_ids or "").replace("|", ";").split(";"):
        pubmed_id = pubmed_id.strip()
        if pubmed_id:
            links.append(f"https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/")
    return ";".join(links)


def number(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def choose_run(report):
    if report.empty:
        return {}

    report = report.copy()
    report["_is_paired"] = report["library_layout"].str.lower().eq("paired")
    report["_fastq_count"] = report["fastq_ftp"].apply(lambda value: len(split_fastqs(value)))
    report["_base_count"] = report["base_count"].apply(number)
    report["_read_count"] = report["read_count"].apply(number)

    # Prefer paired-end reads for downstream assembly, but keep all layouts
    # in all_runs_layouts so this choice can be checked later.
    report = report.sort_values(
        ["_is_paired", "_fastq_count", "_base_count", "_read_count"],
        ascending=False,
    )
    return report.iloc[0].to_dict()


def selection_reason(selected):
    reasons = []
    if str(selected.get("library_layout", "")).lower() == "paired":
        reasons.append("paired_end_preferred")
    else:
        reasons.append("no_paired_run_found")

    if len(split_fastqs(selected.get("fastq_ftp", ""))) > 0:
        reasons.append("has_fastq_links")
    else:
        reasons.append("missing_fastq_links")

    if selected.get("base_count"):
        reasons.append("largest_available_run_tiebreaker")

    return ";".join(reasons)


def all_runs_layouts(report):
    if report.empty:
        return ""
    labels = []
    for _, row in report.iterrows():
        run = row.get("run_accession", "")
        layout = row.get("library_layout", "")
        if run or layout:
            labels.append(f"{run}:{layout}")
    return ";".join(labels)


def all_run_values(report, column):
    if report.empty or column not in report.columns:
        return ""

    values = []
    for _, row in report.iterrows():
        run = row.get("run_accession", "")
        value = row.get(column, "")
        if run or value:
            values.append(f"{run}:{value}")
    return ";".join(values)


def main():
    parser = argparse.ArgumentParser(description="Fetch ENA read metadata.")
    parser.add_argument("--clean-metadata", required=True)
    parser.add_argument("--sequence-metadata", required=True)
    parser.add_argument("--combined-metadata", required=True)
    args = parser.parse_args()

    Path(args.sequence_metadata).parent.mkdir(parents=True, exist_ok=True)
    Path(args.combined_metadata).parent.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(args.clean_metadata)
    rows = []

    for _, sample in metadata.iterrows():
        isolate = sample.get("isolate", "")
        sra_accession = sample.get("SRA_accession", "")

        report = ena_report(sra_accession)
        if report.empty:
            report = query_by_isolate(isolate)

        selected = choose_run(report)
        fastqs = split_fastqs(selected.get("fastq_ftp", ""))

        rows.append(
            {
                "isolate": isolate,
                "BioSample_ID": sample.get("BioSample_ID", ""),
                "assembly_ID": sample.get("assembly_ID", ""),
                "amr_associated_publications": sample.get(
                    "AMR_associated_publications", ""
                ),
                "pubmed_links": pubmed_links(
                    sample.get("AMR_associated_publications", "")
                ),
                "metadata_SRA_accession": sra_accession,
                "BioProject_ID": selected.get("study_accession", ""),
                "BioProject_link": ncbi_link(selected.get("study_accession", "")),
                "SRA_study": selected.get("secondary_study_accession", ""),
                "SRA_study_link": ncbi_link(
                    selected.get("secondary_study_accession", "")
                ),
                "SRA_sample": selected.get("secondary_sample_accession", ""),
                "selected_experiment": selected.get("experiment_accession", ""),
                "selected_run": selected.get("run_accession", ""),
                "selection_reason": selection_reason(selected),
                "all_runs_layouts": all_runs_layouts(report),
                "SRA_submission": selected.get("submission_accession", ""),
                "organism": selected.get("scientific_name", ""),
                "tax_id": selected.get("tax_id", ""),
                "library_layout": selected.get("library_layout", ""),
                "library_name": selected.get("library_name", ""),
                "library_strategy": selected.get("library_strategy", ""),
                "library_source": selected.get("library_source", ""),
                "library_selection": selected.get("library_selection", ""),
                "instrument_platform": selected.get("instrument_platform", ""),
                "instrument_model": selected.get("instrument_model", ""),
                "all_runs_library_layout": all_run_values(report, "library_layout"),
                "all_runs_library_name": all_run_values(report, "library_name"),
                "all_runs_library_strategy": all_run_values(report, "library_strategy"),
                "all_runs_library_source": all_run_values(report, "library_source"),
                "all_runs_library_selection": all_run_values(report, "library_selection"),
                "all_runs_instrument_platform": all_run_values(
                    report, "instrument_platform"
                ),
                "all_runs_instrument_model": all_run_values(report, "instrument_model"),
                "read_count": selected.get("read_count", ""),
                "base_count": selected.get("base_count", ""),
                "read1_fastq_http": to_https(fastqs[0]) if len(fastqs) >= 1 else "",
                "read2_fastq_http": to_https(fastqs[1]) if len(fastqs) >= 2 else "",
                "fastq_ftp": selected.get("fastq_ftp", ""),
                "fastq_md5": selected.get("fastq_md5", ""),
                "submitted_files": selected.get("submitted_ftp", ""),
                "sample_alias": selected.get("sample_alias", ""),
                "sample_description": selected.get("sample_description", ""),
                "center_name": selected.get("center_name", ""),
                "country": selected.get("country", ""),
                "collection_date": selected.get("collection_date", ""),
                "host": selected.get("host", ""),
                "isolation_source": selected.get("isolation_source", ""),
            }
        )

    sequencing_metadata = pd.DataFrame(rows)
    sequencing_metadata.to_csv(args.sequence_metadata, index=False)

    combined = metadata.merge(
        sequencing_metadata,
        on=["BioSample_ID", "assembly_ID", "isolate"],
        how="left",
    )
    combined.to_csv(args.combined_metadata, index=False)


if __name__ == "__main__":
    main()
