#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Create a QQ plot from pyseer output.")
    parser.add_argument("--input", required=True, help="Pyseer output TSV")
    parser.add_argument("--output", required=True, help="Output QQ plot PNG")
    parser.add_argument("--pvalue-column", default="lrt-pvalue")
    args = parser.parse_args()

    table = pd.read_csv(args.input, sep="\t")
    if args.pvalue_column not in table.columns:
        raise SystemExit(f"Missing p-value column: {args.pvalue_column}")

    pvalues = pd.to_numeric(table[args.pvalue_column], errors="coerce")
    pvalues = pvalues[(pvalues > 0) & (pvalues <= 1)].dropna().sort_values()
    if pvalues.empty:
        raise SystemExit("No valid p-values available for QQ plot.")

    observed = -np.log10(pvalues.to_numpy())
    expected = -np.log10(np.arange(1, len(pvalues) + 1) / (len(pvalues) + 1))

    limit = max(float(expected.max()), float(observed.max()))
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(expected, observed, s=12, alpha=0.6, edgecolors="none")
    ax.plot([0, limit], [0, limit], color="black", linewidth=1)
    ax.set_xlabel("Expected -log10(p)")
    ax.set_ylabel("Observed -log10(p)")
    ax.set_title("Pyseer QQ plot")
    ax.set_xlim(0, limit * 1.02)
    ax.set_ylim(0, limit * 1.02)
    fig.tight_layout()
    fig.savefig(args.output, dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
