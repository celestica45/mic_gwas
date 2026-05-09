#!/usr/bin/env python3
import argparse
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def read_unique_patterns(path):
    patterns = set()
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line:
                patterns.add(line)
    return patterns


def format_threshold(value):
    if value is None or pd.isna(value):
        return "NA"
    return f"{value:.3E}"


def summarize_notes(notes):
    counts = Counter()
    for value in notes.fillna("").astype(str):
        value = value.strip()
        if not value:
            counts["no_note"] += 1
            continue
        for token in value.replace(",", ";").split(";"):
            token = token.strip()
            if token:
                counts[token] += 1
    counts.setdefault("bad-chisq", 0)
    counts.setdefault("high-bse", 0)
    counts.setdefault("no_note", 0)
    return counts


def write_threshold_report(path, table, unique_patterns):
    tested = len(table)
    pvalues = table["lrt-pvalue_numeric"].dropna()
    bonferroni = 0.05 / tested if tested else np.nan
    pattern_threshold = 0.05 / len(unique_patterns) if unique_patterns else np.nan
    min_p = pvalues.min() if not pvalues.empty else np.nan
    bonferroni_hits = int((pvalues < bonferroni).sum()) if tested else 0
    pattern_hits = (
        int((pvalues < pattern_threshold).sum()) if unique_patterns else 0
    )
    nominal_hits = int((pvalues < 0.05).sum())

    lines = [
        f"Tested variants: {tested}",
        f"Bonferroni threshold: {format_threshold(bonferroni)}",
        f"Unique patterns: {len(unique_patterns)}",
        f"Pattern threshold: {format_threshold(pattern_threshold)}",
        f"Minimum lrt-pvalue: {format_threshold(min_p)}",
        f"Bonferroni-significant hits: {bonferroni_hits}",
        f"Pattern-significant hits: {pattern_hits}",
        f"Nominal p<0.05 hits: {nominal_hits}",
    ]
    Path(path).write_text("\n".join(lines) + "\n")


def write_notes_summary(path, notes):
    counts = summarize_notes(notes)
    with open(path, "w") as handle:
        handle.write("note\tcount\n")
        for note, count in sorted(counts.items(), key=lambda item: (-item[1], item[0])):
            handle.write(f"{note}\t{count}\n")


def parse_bool(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def shorten_label(value, max_length=35):
    value = str(value)
    if len(value) <= max_length:
        return value
    return value[: max_length - 3] + "..."


def labeled_gene_indices(plot_table, bonferroni, pattern_threshold, label_top_n, label_significant, label_max):
    selected = []

    if label_significant:
        threshold_values = [
            threshold
            for threshold in [bonferroni, pattern_threshold]
            if pd.notna(threshold) and threshold > 0
        ]
        if threshold_values:
            significant_threshold = max(threshold_values)
            significant = plot_table[
                plot_table["lrt-pvalue_numeric"] < significant_threshold
            ].sort_values("lrt-pvalue_numeric")
            selected.extend(significant.index.tolist())

    if label_top_n > 0:
        top_hits = plot_table.nsmallest(label_top_n, "lrt-pvalue_numeric")
        selected.extend(top_hits.index.tolist())

    unique_selected = []
    seen = set()
    for index in selected:
        if index in seen:
            continue
        seen.add(index)
        unique_selected.append(index)

    return unique_selected[: max(label_max, 0)]


def plot_manhattan(
    path,
    table,
    unique_patterns,
    label_top_n,
    label_significant,
    label_max,
):
    plot_table = table.dropna(subset=["lrt-pvalue_numeric"]).copy()
    plot_table = plot_table[
        (plot_table["lrt-pvalue_numeric"] > 0)
        & (plot_table["lrt-pvalue_numeric"] <= 1)
    ]
    if plot_table.empty:
        raise SystemExit("No valid lrt-pvalue values available for Manhattan-like plot.")

    plot_table["minus_log10_p"] = -np.log10(plot_table["lrt-pvalue_numeric"])
    plot_table["index"] = np.arange(1, len(plot_table) + 1)

    tested = len(table)
    bonferroni = 0.05 / tested if tested else np.nan
    pattern_threshold = 0.05 / len(unique_patterns) if unique_patterns else np.nan

    fig, ax = plt.subplots(figsize=(13, 6))
    colors = np.where(plot_table["index"] % 2 == 0, "#2f6f9f", "#d1862f")
    ax.scatter(
        plot_table["index"],
        plot_table["minus_log10_p"],
        s=14,
        c=colors,
        alpha=0.75,
        linewidths=0,
    )

    if pd.notna(bonferroni) and bonferroni > 0:
        ax.axhline(
            -np.log10(bonferroni),
            color="#b22222",
            linestyle="--",
            linewidth=1,
            label="Bonferroni 0.05",
        )
    if pd.notna(pattern_threshold) and pattern_threshold > 0:
        ax.axhline(
            -np.log10(pattern_threshold),
            color="#4b8b3b",
            linestyle=":",
            linewidth=1,
            label="Pattern threshold",
        )

    labels = labeled_gene_indices(
        plot_table,
        bonferroni,
        pattern_threshold,
        label_top_n,
        label_significant,
        label_max,
    )
    for index in labels:
        row = plot_table.loc[index]
        ax.annotate(
            shorten_label(row["variant"]),
            xy=(row["index"], row["minus_log10_p"]),
            xytext=(3, 6),
            textcoords="offset points",
            fontsize=7,
            rotation=30,
            ha="left",
            va="bottom",
        )

    ax.set_xlabel("Gene order in pyseer output")
    ax.set_ylabel("-log10(lrt-pvalue)")
    ax.set_title("Gene GWAS Manhattan-like plot")
    ax.margins(x=0.01)
    ax.legend(loc="best", frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Summarize pyseer gene GWAS thresholds, notes, and plot."
    )
    parser.add_argument("--pyseer", required=True)
    parser.add_argument("--patterns", required=True)
    parser.add_argument("--thresholds", required=True)
    parser.add_argument("--notes", required=True)
    parser.add_argument("--manhattan", required=True)
    parser.add_argument("--label-top-n", type=int, default=10)
    parser.add_argument("--label-significant", default="true")
    parser.add_argument("--label-max", type=int, default=25)
    args = parser.parse_args()

    table = pd.read_csv(args.pyseer, sep="\t")
    if "lrt-pvalue" not in table.columns:
        raise SystemExit("Pyseer output is missing lrt-pvalue column.")
    if "notes" not in table.columns:
        table["notes"] = ""

    table["lrt-pvalue_numeric"] = pd.to_numeric(table["lrt-pvalue"], errors="coerce")
    unique_patterns = read_unique_patterns(args.patterns)

    for output in [args.thresholds, args.notes, args.manhattan]:
        Path(output).parent.mkdir(parents=True, exist_ok=True)

    write_threshold_report(args.thresholds, table, unique_patterns)
    write_notes_summary(args.notes, table["notes"])
    plot_manhattan(
        args.manhattan,
        table,
        unique_patterns,
        max(args.label_top_n, 0),
        parse_bool(args.label_significant),
        max(args.label_max, 0),
    )


if __name__ == "__main__":
    main()
