#!/usr/bin/env python3
"""
Extract summary metrics from RAVEN genome-scale metabolic models (.mat)
AND the number of CDS from annotation .txt files,
and write them to a CSV.

Snakemake mode:
    input.mats: list of .mat files
    input.cds: list of annotation .txt files
    output[0] = output CSV
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
    arr = np.atleast_1d(arr)
    n = 0
    for x in arr:
        if isinstance(x, np.ndarray):
            if x.size == 0:
                continue
            if any(str(y).strip() not in ("", "nan") for y in x.ravel()):
                n += 1
        else:
            s = str(x).strip()
            if s and s.lower() != "nan":
                n += 1
    return n

def _unique_strings_from_mixed(arr) -> set:
    arr = np.atleast_1d(arr)
    unique = set()
    for x in arr:
        if isinstance(x, np.ndarray):
            for y in x.ravel():
                s = str(y).strip()
                if s and s.lower() != "nan":
                    unique.add(s)
        else:
            s = str(x).strip()
            if s and s.lower() != "nan":
                unique.add(s)
    return unique

# ---------- extract CDS count from annotation .txt ----------

def extract_cds_count(txt_path: pathlib.Path) -> int:
    """
    Extract 'number of CDS' from annotation summary file.
    Looks for a line containing something like 'CDS: 2341' or 'Number of CDS = 2341'.
    """
    pattern = re.compile(r"(?:CDS|cds)[^\d]*(\d+)")
    with open(txt_path, "r", encoding="utf-8") as f:
        for line in f:
            match = pattern.search(line)
            if match:
                return int(match.group(1))
    return np.nan  # in case nothing matches

# ---------- exchange / transport detection ----------

def _detect_external_comps(model) -> set:
    if not hasattr(model, "comps") or not hasattr(model, "compNames"):
        return set()

    comp_ids = np.atleast_1d(model.comps)
    comp_names = np.atleast_1d(model.compNames)

    ext = []
    for i, (cid, cname) in enumerate(zip(comp_ids, comp_names)):
        cid = str(cid).lower()
        cname = str(cname).lower()

        if any(term in cname for term in ["extracellular", "external", "periplasm", "outside"]):
            ext.append(i + 1)
        elif cid in ("e", "e0", "extracellular"):
            ext.append(i + 1)

    return set(ext)

# ---------- core metrics computation ----------

def compute_metrics(model, mat_path: pathlib.Path):
    stats = {}

    # file-level identifiers
    stats["file"] = str(mat_path)
    stats["strain"] = mat_path.stem.split("_")[0]   # ESL0001_draft.mat → ESL0001

    # GEM structural info
    stats["n_rxns"] = int(getattr(model, "rxns").size)
    stats["n_mets"] = int(getattr(model, "mets").size)
    stats["n_genes"] = int(getattr(model, "genes", np.array([])).size)

    # reversibility
    if hasattr(model, "rev"):
        rev = np.asarray(model.rev).ravel().astype(int)
        stats["n_rev"] = int((rev == 1).sum())
        stats["n_irrev"] = int(rev.size - stats["n_rev"])
        stats["frac_rev"] = stats["n_rev"] / rev.size
    else:
        stats["n_rev"] = stats["n_irrev"] = stats["frac_rev"] = np.nan

    # gene rules
    if hasattr(model, "grRules"):
        n_with = _non_empty_reaction_annotations(model.grRules)
        stats["n_rxns_with_gr"] = n_with
        stats["frac_rxns_with_gr"] = n_with / stats["n_rxns"]
    else:
        stats["n_rxns_with_gr"] = stats["frac_rxns_with_gr"] = np.nan

    # EC numbers
    if hasattr(model, "eccodes"):
        n_with = _non_empty_reaction_annotations(model.eccodes)
        stats["n_rxns_with_EC"] = n_with
        stats["frac_rxns_with_EC"] = n_with / stats["n_rxns"]
    else:
        stats["n_rxns_with_EC"] = stats["frac_rxns_with_EC"] = np.nan

    # subsystems
    if hasattr(model, "subSystems"):
        n_with = _non_empty_reaction_annotations(model.subSystems)
        stats["n_rxns_with_subsystem"] = n_with
        stats["frac_rxns_with_subsystem"] = n_with / stats["n_rxns"]
        stats["n_unique_subsystems"] = len(_unique_strings_from_mixed(model.subSystems))
    else:
        stats["n_rxns_with_subsystem"] = stats["frac_rxns_with_subsystem"] = np.nan
        stats["n_unique_subsystems"] = np.nan

    # exchange/transport
    S = sp.csc_matrix(getattr(model, "S"))
    metComps = getattr(model, "metComps", None)
    ext = _detect_external_comps(model)

    n_exch = 0
    n_trans = 0

    for j in range(stats["n_rxns"]):
        rows = S.indices[S.indptr[j]:S.indptr[j+1]]
        if rows.size == 0:
            continue

        # transport
        if metComps is not None:
            comps = np.atleast_1d(metComps)[rows]
            if len(set(int(c) for c in comps)) >= 2:
                n_trans += 1

        # exchange
        if rows.size == 1 and metComps is not None:
            comp = int(np.atleast_1d(metComps)[rows[0]])
            if comp in ext:
                n_exch += 1

    stats["n_exchange_rxns"] = n_exch
    stats["frac_exchange_rxns"] = n_exch / stats["n_rxns"]
    stats["n_transport_rxns"] = n_trans
    stats["frac_transport_rxns"] = n_trans / stats["n_rxns"]

    return stats

# ---------- Snakemake / CLI wrappers ----------

def load_raven_model(mat_path: pathlib.Path):
    d = loadmat(mat_path, squeeze_me=True, struct_as_record=False)
    key = next(k for k in d.keys() if not k.startswith("__"))
    return d[key]

def run(mat_paths, cds_paths, out_csv: pathlib.Path):
    rows = []
    for mat, cds in zip(mat_paths, cds_paths):
        model = load_raven_model(mat)
        stats = compute_metrics(model, mat)

        # Add CDS count
        stats["n_cds"] = extract_cds_count(cds)

        rows.append(stats)

    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)

# ---------- CLI entry ----------

def _parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mat", nargs="+", required=True, help="MAT model files")
    ap.add_argument("--cds", nargs="+", required=True, help="Annotation TXT files")
    ap.add_argument("-o", "--output", required=True)
    return ap.parse_args()

if __name__ == "__main__":
    if "snakemake" in globals():
        mats = [pathlib.Path(p) for p in snakemake.input.mats]
        cds = [pathlib.Path(p) for p in snakemake.input.cds]
        out = pathlib.Path(snakemake.output[0])
        run(mats, cds, out)

    else:
        args = _parse_args()
        mats = [pathlib.Path(p) for p in args.mat]
        cds = [pathlib.Path(p) for p in args.cds]
        run(mats, cds, pathlib.Path(args.output))
