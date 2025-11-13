#!/usr/bin/env python3
"""
Combine plate-reader Excel files + one metadata.xlsx into a single TSV.

Output columns:
  species, replicate, media, cond, time, value
"""

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, Optional, List

import pandas as pd
import numpy as np

# ---------- Config / heuristics ----------
COND_PATTERNS = [
    (re.compile(r"microaer", re.I), "microaerobic"),
    (re.compile(r"anaer", re.I), "anaerobic"),
]
WELL_REGEX_96 = re.compile(r"^[A-H](?:[1-9]|1[0-2])$")

def find_condition_from_filename(path: Path) -> Optional[str]:
    name = path.stem.lower()
    for pat, label in COND_PATTERNS:
        if pat.search(name):
            return label
    return None

def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(c).strip().lower() for c in df.columns]
    return df

def to_minutes(x) -> Optional[float]:
    if pd.isna(x):
        return None
    if isinstance(x, (int, float, np.integer, np.floating)):
        return float(x)
    if isinstance(x, (pd.Timedelta, np.timedelta64)):
        return float(pd.to_timedelta(x).total_seconds()) / 60.0
    s = str(x).strip()
    if re.fullmatch(r"\d{1,2}:\d{2}:\d{2}", s):
        h, m, sec = [int(p) for p in s.split(":")]
        return h * 60 + m + sec / 60.0
    if re.fullmatch(r"\d{1,2}:\d{2}", s):
        m, sec = [int(p) for p in s.split(":")]
        return m + sec / 60.0
    m = re.fullmatch(r"\s*([0-9]*\.?[0-9]+)\s*([a-zA-Z]+)\s*", s)
    if m:
        val = float(m.group(1)); unit = m.group(2).lower()
        if unit.startswith("s"): return val / 60.0
        if unit.startswith("min") or unit == "m": return val
        if unit.startswith("h"): return val * 60.0
    try:
        return float(s)
    except Exception:
        return None

# ---------- Metadata ----------
def read_metadata(metadata_path: Path) -> Dict[str, pd.DataFrame]:
    xl = pd.ExcelFile(metadata_path)
    out = {}
    for sheet in xl.sheet_names:
        df = xl.parse(sheet)
        df = normalize_columns(df)
        if "species" not in df.columns and "bacteria" in df.columns:
            df = df.rename(columns={"bacteria": "species"})
        if "well" in df.columns:
            df["well"] = (
                df["well"].astype(str).str.strip().str.upper()
                .str.replace(r"^([A-H])0?([1-9]|1[0-2])$", r"\1\2", regex=True)
            )
        keep = [c for c in ["well", "species", "media", "replicate"] if c in df.columns]
        out[sheet] = df[keep].copy()
    return out

def best_metadata_sheet_for_file(meta_by_sheet: Dict[str, pd.DataFrame], file_path: Path) -> Optional[str]:
    fname = file_path.stem.lower()
    for s in meta_by_sheet:
        if s.lower() == fname:
            return s
    candidates = [s for s in meta_by_sheet if s.lower() in fname or fname in s.lower()]
    if len(candidates) == 1:
        return candidates[0]
    if candidates:
        def score(s):
            sl = s.lower()
            return sum(a == b for a, b in zip(sl, fname))
        return sorted(candidates, key=score, reverse=True)[0]
    return list(meta_by_sheet.keys())[0] if meta_by_sheet else None

# ---------- Plate parsing (header-scan) ----------
def parse_plate_file(path: Path) -> pd.DataFrame:
    """
    Detect the header row that contains 'Time' and at least one well (A1..H12).
    Keep 'Time' + well columns. Melt to (well, time, value).
    """
    xl = pd.ExcelFile(path)
    last_err = None
    for sheet in xl.sheet_names:
        try:
            raw = xl.parse(sheet, header=None)
            header_row = None
            for i in range(len(raw)):
                row = raw.iloc[i].astype(str).str.strip().tolist()
                has_time = any(cell.lower() == "time" for cell in row)
                has_well = any(WELL_REGEX_96.fullmatch(cell) for cell in row)
                if has_time and has_well:
                    header_row = i
                    break
            if header_row is None:
                continue

            header = raw.iloc[header_row].astype(str).str.strip().tolist()
            data = raw.iloc[header_row + 1:].copy()
            if len(header) < data.shape[1]:
                header = header + [f"extra_{j}" for j in range(data.shape[1] - len(header))]
            data.columns = header[: data.shape[1]]

            time_col = next((c for c in data.columns if str(c).lower() == "time"), None)
            if not time_col:
                raise ValueError("No 'Time' column after header detection.")
            well_cols = [c for c in data.columns if WELL_REGEX_96.fullmatch(str(c))]
            if not well_cols:
                raise ValueError("No well columns (A1..H12) found.")

            df = data[[time_col] + well_cols].copy()
            df = df.dropna(how="all")
            df["time"] = df[time_col].apply(to_minutes)
            df = df.dropna(subset=["time"])

            long = df.melt(id_vars=["time"], value_vars=well_cols, var_name="well", value_name="value")
            long["value"] = pd.to_numeric(long["value"], errors="coerce")
            long = long.dropna(subset=["value"]).reset_index(drop=True)
            long["_sheet"] = sheet
            return long
        except Exception as e:
            last_err = e
            continue
    raise ValueError(f"{path.name}: Could not parse any sheet. Last error: {last_err}")

# ---------- Replicates across ALL files ----------
def assign_replicates(big: pd.DataFrame) -> pd.DataFrame:
    """
    If metadata provides 'replicate', normalize to 1..n per (species, media, cond).
    Otherwise derive replicate IDs by order-of-first-appearance of 'well' per group.
    """
    if "replicate" in big.columns and big["replicate"].notna().any():
        big["replicate"] = (
            big.groupby(["species", "media", "cond"], group_keys=False)["replicate"]
               .transform(lambda s: pd.factorize(s, sort=False)[0] + 1)
               .astype(int)
        )
        return big

    # First index when each (species, media, cond, file, well) appears (so wells from
    # different files become distinct replicates if they’re the same physical well id)
    tmp = big.reset_index()
    firsts = (
        tmp.groupby(["species", "media", "cond", "_file", "well"], as_index=False)
           .agg(first_idx=("index", "min"))
    )
    # Within each (species, media, cond), rank by first appearance order:
    firsts["replicate"] = (
        firsts.sort_values(["species","media","cond","first_idx"])
              .groupby(["species","media","cond"])
              .cumcount() + 1
    )
    big = big.merge(
        firsts.drop(columns=["first_idx"]),
        on=["species", "media", "cond", "_file", "well"],
        how="left",
        validate="many_to_one",
    )
    return big

# ---------- Main ----------
def main():
    ap = argparse.ArgumentParser()
    # Accept many --raw flags (Snakemake will pass a list)
    ap.add_argument("--raw", action="append", required=True,
                    help="Path to one experiment .xlsx file. Use multiple --raw to combine many files.")
    ap.add_argument("--metadata", required=True, help="Path to metadata.xlsx")
    ap.add_argument("--output", required=True, help="Output TSV path (single combined file)")
    args = ap.parse_args()

    raw_paths = [Path(p) for p in args.raw]
    metadata_path = Path(args.metadata)
    out_path = Path(args.output)

    missing = [str(p) for p in raw_paths if not p.exists()]
    if missing:
        print("ERROR: missing input files:\n  " + "\n  ".join(missing), file=sys.stderr)
        sys.exit(1)
    if not metadata_path.exists():
        print(f"ERROR: metadata file not found: {metadata_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Reading metadata: {metadata_path.name}")
    meta_by_sheet = read_metadata(metadata_path)
    if not meta_by_sheet:
        print("ERROR: No sheets found in metadata.xlsx.", file=sys.stderr)
        sys.exit(1)

    rows = []
    errors = []
    for xlsx in sorted(raw_paths):
        try:
            cond = find_condition_from_filename(xlsx) or "unknown"
            plate_long = parse_plate_file(xlsx)  # well, time, value, _sheet
            meta_sheet = best_metadata_sheet_for_file(meta_by_sheet, xlsx)
            meta_df = meta_by_sheet.get(meta_sheet)
            if meta_df is None or meta_df.empty:
                raise ValueError(f"No usable metadata sheet found (sheet='{meta_sheet}')")

            for col in ("well", "species", "media"):
                if col not in meta_df.columns:
                    raise ValueError(f"Metadata sheet '{meta_sheet}' must contain '{col}' column.")

            merged = plate_long.merge(
                meta_df,
                on="well",
                how="inner",
                validate="many_to_one"
            )
            merged["_file"] = xlsx.name
            merged["cond"] = cond
            merged = merged[pd.to_numeric(merged["value"], errors="coerce").notna()].copy()
            merged["value"] = merged["value"].astype(float)
            rows.append(merged)
            print(f"Parsed {xlsx.name} (wells merged: {merged['well'].nunique()}, rows: {len(merged)})")
        except Exception as e:
            errors.append(f"{xlsx.name}: {e}")
            print(f"WARNING: {xlsx.name}: {e}", file=sys.stderr)

    if not rows:
        print("ERROR: No data parsed from experiment files.", file=sys.stderr)
        if errors:
            print("\nDetails:", file=sys.stderr)
            for err in errors:
                print(" - " + err, file=sys.stderr)
        sys.exit(1)

    big = pd.concat(rows, ignore_index=True)
    big = assign_replicates(big)

    # --- normalize time values ---
    # Round to nearest 0.5 min, then force uniform spacing per file if needed
    big["time"] = big["time"].round(1)

    # Reindex times per file if they look irregular
    def regularize_time(df):
        times = np.sort(df["time"].unique())
        # Check if spacing is roughly constant
        if len(times) > 1:
            diffs = np.diff(times)
            step = np.median(diffs)
            if step < 0.1:
                step = 30.0  # fallback default 30 min if near-zero
            # force 0-based evenly spaced sequence
            new_times = np.arange(0, step * len(times), step)
            mapping = dict(zip(times, new_times))
            df["time"] = df["time"].map(mapping)
        else:
            df["time"] = 0.0
        return df

    big = big.groupby("_file", group_keys=False).apply(regularize_time)

    # --- prepare final output ---
    final = big[["species", "replicate", "media", "cond", "time", "value"]].copy()
    # convert to int if all times are whole numbers
    if np.allclose(final["time"], final["time"].round()):
        final["time"] = final["time"].astype(int)
    final = final.sort_values(["species", "media", "cond", "replicate", "time"]).reset_index(drop=True)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(out_path, sep="\t", index=False)
    print(f"\nWrote: {out_path}  (rows: {len(final)}, species: {final['species'].nunique()})")

if __name__ == "__main__":
    pd.options.mode.copy_on_write = True
    main()
