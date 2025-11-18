#!/usr/bin/env python3

"""
Summarize GEM reconstruction metrics from a metrics CSV
(e.g. draft_metrics.csv or relaxed_metrics.csv).

Usage:
    python model_statistics.py draft_metrics.csv
    python model_statistics.py relaxed_metrics.csv

Required columns in the CSV:
    - strain   : identifier for the strain / model
    - n_rxns   : number of reactions in the model
    - n_mets   : number of metabolites
    - n_genes  : number of genes present in the model (metabolic genes)
    - n_cds    : total number of CDS predicted in the genome

The script computes, for each model:
    - frac_metabolic_genes      = n_genes / n_cds
    - metabolites_per_reaction  = n_mets / n_rxns
    - genes_per_reaction_global = n_genes / n_rxns

and prints:
    1) A per-model table
    2) A global summary (mean, median, std, min, max) across all models
"""

import argparse
import pathlib
import pandas as pd
import numpy as np


def parse_args():
    ap = argparse.ArgumentParser(
        description="Summarize GEM metrics from a CSV."
    )
    ap.add_argument(
        "metrics_csv",
        help="Path to metrics CSV (e.g. draft_metrics.csv).",
    )
    ap.add_argument(
        "-o", "--output",
        help="Output Markdown file (default: <csvname>_report.md)",
        default=None
    )
    return ap.parse_args()


def main():
    args = parse_args()
    metrics_path = pathlib.Path(args.metrics_csv)

    # Determine output markdown path
    if args.output is None:
        md_path = metrics_path.with_name(metrics_path.stem + "_report.md")
    else:
        md_path = pathlib.Path(args.output)

    df = pd.read_csv(metrics_path)

    # Required columns
    required_cols = ["strain", "n_rxns", "n_mets", "n_genes", "n_cds"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in CSV: {missing}")

    # Compute derived metrics
    df["frac_metabolic_genes"] = df["n_genes"] / df["n_cds"]
    df["metabolites_per_reaction"] = df["n_mets"] / df["n_rxns"]
    df["genes_per_reaction_global"] = df["n_genes"] / df["n_rxns"]

    # Per-model metrics
    per_model_cols = [
        "strain",
        "n_rxns",
        "n_mets",
        "n_genes",
        "n_cds",
        "frac_metabolic_genes",
        "frac_rev",
        "frac_rxns_with_gr",
        "frac_rxns_with_EC",
        "frac_rxns_with_subsystem",
        "metabolites_per_reaction",
        "genes_per_reaction_global",
        "n_exchange_rxns",
        "n_transport_rxns",
    ]
    per_model_cols = [c for c in per_model_cols if c in df.columns]
    per_model = df[per_model_cols].copy()

    # Global summary (mean, median, etc.)
    metrics_for_summary = [
        "n_rxns",
        "n_mets",
        "n_genes",
        "n_cds",
        "frac_metabolic_genes",
        "frac_rev",
        "frac_rxns_with_gr",
        "frac_rxns_with_EC",
        "frac_rxns_with_subsystem",
        "metabolites_per_reaction",
        "genes_per_reaction_global",
        "n_exchange_rxns",
        "n_transport_rxns",
    ]
    metrics_for_summary = [c for c in metrics_for_summary if c in df.columns]

    summary = {}
    for col in metrics_for_summary:
        series = df[col].dropna()
        if series.empty:
            continue
        summary[col] = {
            "mean": float(series.mean()),
            "median": float(series.median()),
            "std": float(series.std()),
            "min": float(series.min()),
            "max": float(series.max()),
        }

    summary_df = pd.DataFrame(summary).T[
        ["mean", "median", "std", "min", "max"]
    ]

    # Print to console
    print("\n=== Per-model summary ===\n")
    print(per_model.to_string(index=False))

    print("\n=== Global summary across all models ===\n")
    print(summary_df.to_string())

    # ---- Write Markdown report ----
    with open(md_path, "w", encoding="utf-8") as md:
        md.write(f"# GEM Reconstruction Statistics Report\n")
        md.write(f"Source CSV: `{metrics_path.name}`\n\n")

        # Per-model table
        md.write("## Per-model summary\n\n")
        md.write(per_model.to_markdown(index=False))
        md.write("\n\n")

        # Global summary
        md.write("## Global statistics\n\n")
        md.write(summary_df.to_markdown())
        md.write("\n\n")

    print(f"\nMarkdown report written to: {md_path}\n")


if __name__ == "__main__":
    main()
