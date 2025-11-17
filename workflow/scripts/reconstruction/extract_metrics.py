#!/usr/bin/env python3
"""
Extract summary metrics from RAVEN genome-scale metabolic models (.mat)
and write them to a CSV.

Works both as:
  - a standalone script:
        python extract_gem_metrics.py models/*.mat -o gem_metrics.csv
  - a Snakemake `script:`:
        script: "workflow/scripts/extract_gem_metrics.py"
    where snakemake.input = list of .mat files
          snakemake.output[0] = output CSV path
"""

import argparse
import pathlib
import re

import numpy as np
import pandas as pd
from scipy.io import loadmat
import scipy.sparse as sp


# ---------- helpers for RAVEN-style mixed cell arrays ----------

def _non_empty_reaction_annotations(arr) -> int:
    """
    Count entries that have at least one non-empty string.
    Handles RAVEN-style arrays where entries can be:
      - string
      - numpy array of strings
      - empty arrays
    """
    arr = np.atleast_1d(arr)
    n = 0
    for x in arr:
        if isinstance(x, np.ndarray):
            if x.size == 0:
                continue
            has_any = False
            for y in x.ravel():
                s = str(y).strip()
                if s and s.lower() not in ("nan",):
                    has_any = True
                    break
            if has_any:
                n += 1
        else:
            s = str(x).strip()
            if s and s.lower() not in ("nan",):
                n += 1
    return n


def _unique_strings_from_mixed(arr) -> set:
    """Collect unique non-empty strings from mixed RAVEN-style arrays."""
    arr = np.atleast_1d(arr)
    unique = set()
    for x in arr:
        if isinstance(x, np.ndarray):
            for y in x.ravel():
                s = str(y).strip()
                if s and s.lower() not in ("nan",):
                    unique.add(s)
        else:
            s = str(x).strip()
            if s and s.lower() not in ("nan",):
                unique.add(s)
    return unique


# ---------- exchange / transport detection ----------

def _detect_external_comps(model) -> set:
    """
    Heuristic: identify external compartments by ID or name.
    Returns a set of compartment indices (1-based, as in RAVEN's metComps).
    """
    if not hasattr(model, "comps") or not hasattr(model, "compNames"):
        return set()

    comp_ids = np.atleast_1d(model.comps)
    comp_names = np.atleast_1d(model.compNames)

    external_idxs = []
    for i, (cid, cname) in enumerate(zip(comp_ids, comp_names)):
        cid_str = str(cid).lower()
        cname_str = str(cname).lower()

        # Name-based
        if any(word in cname_str for word in
               ["extracellular", "external", "periplasm", "outside", "boundary"]):
            external_idxs.append(i + 1)  # MATLAB 1-based indices
        # ID-based
        elif cid_str in ("e", "e0", "extracellular"):
            external_idxs.append(i + 1)

    return set(external_idxs)


# ---------- core metrics computation ----------

def compute_metrics(model, mat_path: pathlib.Path):
    stats = {}

    # file-level identifiers
    stats["file"] = str(mat_path)
    stats["id"] = str(getattr(model, "id", ""))
    stats["name"] = str(getattr(model, "name", ""))

    # basic sizes
    stats["n_rxns"] = int(getattr(model, "rxns").size)
    stats["n_mets"] = int(getattr(model, "mets").size)
    stats["n_genes"] = int(getattr(model, "genes", np.array([])).size)

    # reversibility
    if hasattr(model, "rev"):
        rev = np.asarray(model.rev).ravel().astype(int)
        stats["n_rev"] = int(np.sum(rev == 1))
        stats["n_irrev"] = int(rev.size - stats["n_rev"])
        stats["frac_rev"] = stats["n_rev"] / rev.size if rev.size > 0 else np.nan
    else:
        stats["n_rev"] = stats["n_irrev"] = stats["frac_rev"] = np.nan

    # gene rules (grRules)
    if hasattr(model, "grRules"):
        n_with_gr = _non_empty_reaction_annotations(model.grRules)
        stats["n_rxns_with_gr"] = int(n_with_gr)
        stats["frac_rxns_with_gr"] = (
            n_with_gr / stats["n_rxns"] if stats["n_rxns"] > 0 else np.nan
        )
    else:
        stats["n_rxns_with_gr"] = stats["frac_rxns_with_gr"] = np.nan

    # genes per reaction from rxnGeneMat
    if hasattr(model, "rxnGeneMat"):
        rg = model.rxnGeneMat
        # scipy loads MATLAB sparse matrices as scipy.sparse.spmatrix
        gene_counts = np.array((rg != 0).sum(axis=1)).ravel()
        if gene_counts.size > 0:
            stats["mean_genes_per_rxn"] = float(gene_counts.mean())
            stats["median_genes_per_rxn"] = float(np.median(gene_counts))
        else:
            stats["mean_genes_per_rxn"] = stats["median_genes_per_rxn"] = np.nan
    else:
        stats["mean_genes_per_rxn"] = stats["median_genes_per_rxn"] = np.nan

    # subsystems
    if hasattr(model, "subSystems"):
        n_with_ss = _non_empty_reaction_annotations(model.subSystems)
        unique_ss = _unique_strings_from_mixed(model.subSystems)
        stats["n_rxns_with_subsystem"] = int(n_with_ss)
        stats["frac_rxns_with_subsystem"] = (
            n_with_ss / stats["n_rxns"] if stats["n_rxns"] > 0 else np.nan
        )
        stats["n_unique_subsystems"] = len(unique_ss)
    else:
        stats["n_rxns_with_subsystem"] = stats["frac_rxns_with_subsystem"] = np.nan
        stats["n_unique_subsystems"] = np.nan

    # EC numbers
    if hasattr(model, "eccodes"):
        n_with_EC = _non_empty_reaction_annotations(model.eccodes)
        unique_ec = _unique_strings_from_mixed(model.eccodes)
        stats["n_rxns_with_EC"] = int(n_with_EC)
        stats["frac_rxns_with_EC"] = (
            n_with_EC / stats["n_rxns"] if stats["n_rxns"] > 0 else np.nan
        )
        stats["n_unique_EC"] = len(unique_ec)
    else:
        stats["n_rxns_with_EC"] = stats["frac_rxns_with_EC"] = np.nan
        stats["n_unique_EC"] = np.nan

    # exchange / transport using S and metComps
    S = getattr(model, "S")
    if not isinstance(S, sp.spmatrix):
        S = sp.csc_matrix(S)
    else:
        S = S.tocsc()

    metComps = getattr(model, "metComps", None)
    n_rxns = stats["n_rxns"]

    external_comps = set()
    if hasattr(model, "comps") and hasattr(model, "compNames") and hasattr(model, "metComps"):
        external_comps = _detect_external_comps(model)

    n_exch = 0
    n_trans = 0

    for j in range(n_rxns):
        start, end = S.indptr[j], S.indptr[j + 1]
        rows = S.indices[start:end]
        if rows.size == 0:
            continue

        # transport: metabolites from >= 2 different compartments
        if metComps is not None:
            comps = np.atleast_1d(metComps)[rows]
            uniq_comps = set(int(c) for c in comps)
            if len(uniq_comps) >= 2:
                n_trans += 1

        # exchange: single metabolite in an "external" compartment (heuristic)
        if external_comps and metComps is not None and rows.size == 1:
            comp = int(np.atleast_1d(metComps)[rows[0]])
            if comp in external_comps:
                n_exch += 1

    stats["n_exchange_rxns"] = int(n_exch)
    stats["frac_exchange_rxns"] = n_exch / n_rxns if n_rxns > 0 else np.nan
    stats["n_transport_rxns"] = int(n_trans)
    stats["frac_transport_rxns"] = n_trans / n_rxns if n_rxns > 0 else np.nan

    return stats


# ---------- I/O wrappers ----------

def load_raven_model(mat_path: pathlib.Path):
    """
    Load a .mat file and return the first non-internal variable
    (assumed to be the RAVEN model struct).
    """
    data = loadmat(mat_path, squeeze_me=True, struct_as_record=False)
    model_key = next(k for k in data.keys() if not k.startswith("__"))
    return data[model_key]


def run(mat_paths, out_csv: pathlib.Path):
    rows = []
    for p in mat_paths:
        model = load_raven_model(p)
        rows.append(compute_metrics(model, p))

    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)


# ---------- CLI entry point ----------

def _parse_args():
    ap = argparse.ArgumentParser(
        description="Extract metrics from RAVEN GEM .mat models into a CSV."
    )
    ap.add_argument(
        "input",
        nargs="+",
        help="One or more .mat files (wildcards allowed if your shell supports them).",
    )
    ap.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output CSV file.",
    )
    return ap.parse_args()


if __name__ == "__main__":
    # Snakemake mode
    if "snakemake" in globals():
        mat_paths = [pathlib.Path(p) for p in snakemake.input]  # type: ignore
        out_csv = pathlib.Path(snakemake.output[0])             # type: ignore
        run(mat_paths, out_csv)

    # Normal CLI mode
    else:
        args = _parse_args()
        mat_paths = [pathlib.Path(p) for p in args.input]
        out_csv = pathlib.Path(args.output)
        run(mat_paths, out_csv)
