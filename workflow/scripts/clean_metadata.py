#!/usr/bin/env python3
import argparse
import math
import re
from pathlib import Path

import pandas as pd


def has_parseable_mic(value):
    return parse_mic_value(value) is not None


def parse_mic_value(value):
    if pd.isna(value):
        return None
    match = re.search(r"\d+(\.\d+)?", str(value))
    if match is None:
        return None
    mic = float(match.group(0))
    if mic <= 0:
        return None
    return mic


def write_mic_log2_outputs(cleaned, stats_output, counts_output, pie_chart_output):
    Path(stats_output).parent.mkdir(parents=True, exist_ok=True)
    Path(counts_output).parent.mkdir(parents=True, exist_ok=True)
    Path(pie_chart_output).parent.mkdir(parents=True, exist_ok=True)

    mic_log2 = cleaned["mic_log2"].dropna()
    stats = mic_log2.describe(percentiles=[0.25, 0.5, 0.75]).rename(
        {
            "count": "n",
            "25%": "q1",
            "50%": "median",
            "75%": "q3",
        }
    )
    stats.to_frame(name="mic_log2").reset_index(names="statistic").to_csv(
        stats_output, index=False
    )

    counts = (
        cleaned["mic_log2"]
        .value_counts()
        .sort_index()
        .rename_axis("mic_log2")
        .reset_index(name="count")
    )
    counts["label"] = counts.apply(
        lambda row: f"log2 MIC {row['mic_log2']:g} (n={int(row['count'])})",
        axis=1,
    )
    counts.to_csv(counts_output, index=False)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    total = counts["count"].sum()
    rare_threshold = 20 if total >= 100 else 1
    common = counts[counts["count"] >= rare_threshold].copy()
    rare = counts[counts["count"] < rare_threshold]
    if not rare.empty:
        common = pd.concat(
            [
                common,
                pd.DataFrame(
                    {
                        "mic_log2": ["Other"],
                        "count": [rare["count"].sum()],
                    }
                ),
            ],
            ignore_index=True,
        )
    common["plot_label"] = common.apply(
        lambda row: f"{row['mic_log2']:g}: {int(row['count'])}"
        if row["mic_log2"] != "Other"
        else f"Other: {int(row['count'])}",
        axis=1,
    )

    colors = [
        "#4E79A7",
        "#F28E2B",
        "#59A14F",
        "#E15759",
        "#76B7B2",
        "#B07AA1",
        "#9C755F",
    ]
    fig, ax = plt.subplots(figsize=(10, 7))
    wedges, _, autotexts = ax.pie(
        common["count"],
        labels=None,
        colors=colors[: len(common)],
        autopct=lambda pct: f"{pct:.1f}%" if pct >= 4 else "",
        startangle=90,
        counterclock=False,
        pctdistance=0.72,
        wedgeprops={"linewidth": 1, "edgecolor": "white"},
        textprops={"fontsize": 10, "color": "white", "weight": "bold"},
    )
    ax.legend(
        wedges,
        common["plot_label"],
        title="MIC log2: count",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
    )
    ax.set_title(f"MIC log2 distribution (n={int(total)})")
    ax.axis("equal")
    for autotext in autotexts:
        autotext.set_fontsize(10)
    fig.tight_layout(rect=[0, 0, 0.78, 1])
    fig.savefig(pie_chart_output, dpi=200)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Clean raw AMR metadata.")
    parser.add_argument("--input", required=True, help="Raw metadata CSV")
    parser.add_argument("--output", required=True, help="Clean metadata CSV")
    parser.add_argument(
        "--output-with-assembly",
        required=True,
        help="Clean metadata CSV keeping only rows with assembly_ID",
    )
    parser.add_argument("--summary", required=True, help="Metadata cleaning summary")
    parser.add_argument("--mic-log2-stats", required=True, help="MIC log2 statistics CSV")
    parser.add_argument("--mic-log2-counts", required=True, help="MIC log2 counts CSV")
    parser.add_argument("--mic-log2-pie-chart", required=True, help="MIC log2 pie chart PNG")
    parser.add_argument("--antibiotic-name", required=True)
    parser.add_argument("--organism", required=True)
    args = parser.parse_args()

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output_with_assembly).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(args.input)

    # First clean column names so phenotype-gen_measurement becomes gen_measurement.
    metadata = metadata.rename(
        columns=lambda column: column.removeprefix("phenotype-")
    )

    if "gen_measurement" not in metadata.columns:
        raise SystemExit("Missing required column after cleaning: gen_measurement")
    if "assembly_ID" not in metadata.columns:
        raise SystemExit("Missing required column after cleaning: assembly_ID")
    if "antibiotic_name" not in metadata.columns:
        raise SystemExit("Missing required column after cleaning: antibiotic_name")
    if "organism" not in metadata.columns:
        raise SystemExit("Missing required column after cleaning: organism")

    # Start with only the antibiotic and organism this dataset is for.
    filtered = metadata[
        (metadata["antibiotic_name"].astype(str).str.lower() == args.antibiotic_name.lower())
        & (metadata["organism"].astype(str) == args.organism)
    ].copy()

    rows_after_filter = len(filtered)
    filtered["mic_value"] = filtered["gen_measurement"].apply(parse_mic_value)
    missing_or_unparseable_mic = filtered["mic_value"].isna()

    drop_mask = missing_or_unparseable_mic
    cleaned = filtered[~drop_mask].copy()
    cleaned["mic_log2"] = cleaned["mic_value"].apply(math.log2)

    # Keep the antibiotic name from config and place it first.
    cleaned["antibiotic_name"] = args.antibiotic_name
    cleaned = cleaned[
        ["antibiotic_name"]
        + [column for column in cleaned.columns if column != "antibiotic_name"]
    ]

    has_assembly = cleaned["assembly_ID"].notna() & (
        cleaned["assembly_ID"].astype(str).str.strip() != ""
    )
    cleaned_with_assembly = cleaned[has_assembly].copy()
    cleaned.to_csv(args.output, index=False)
    cleaned_with_assembly.to_csv(args.output_with_assembly, index=False)
    write_mic_log2_outputs(
        cleaned_with_assembly,
        args.mic_log2_stats,
        args.mic_log2_counts,
        args.mic_log2_pie_chart,
    )

    with open(args.summary, "w") as handle:
        handle.write(f"antibiotic: {args.antibiotic_name}\n")
        handle.write(f"organism: {args.organism}\n")
        handle.write(f"rows_after_antibiotic_organism_filter: {rows_after_filter}\n")
        handle.write(f"rows_kept: {len(cleaned)}\n")
        handle.write(f"rows_dropped: {int(drop_mask.sum())}\n")
        handle.write(f"rows_with_assembly: {len(cleaned_with_assembly)}\n")
        handle.write(f"rows_without_assembly: {int((~has_assembly).sum())}\n")
        handle.write("\n")
        handle.write("drop_reasons:\n")
        handle.write(
            f"  missing_or_unparseable_mic: {int(missing_or_unparseable_mic.sum())}\n"
        )


if __name__ == "__main__":
    main()
