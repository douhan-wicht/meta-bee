#!/usr/bin/env python3

"""
growth_amiga.py
-----------------------
Fit AMiGA Growth Models to plate-reader data.

Input TSV must contain:
    species, media, cond, replicate, time, value

Output: CSV table with
    species, media, cond, replicate,
    mumax, lag, max_derivative, growth_call
"""

import argparse
import pandas as pd
import numpy as np
import sys
import os
import amiga

# ----------------------------------------------
# Parse arguments
# ----------------------------------------------
parser = argparse.ArgumentParser(description="AMiGA Growth Curve Fitting")
parser.add_argument("--input-tsv", required=True)
parser.add_argument("--output-csv", required=True)
args = parser.parse_args()

input_tsv = args.input_tsv
output_csv = args.output_csv

print("Input :", input_tsv)
print("Output:", output_csv)


# ----------------------------------------------
# Load data
# ----------------------------------------------
df = pd.read_csv(input_tsv, sep="\t")

required = {"species", "media", "cond", "replicate", "time", "value"}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns: {missing}")

df = df[df["value"] > 0].copy()
df = df[np.isfinite(df["time"]) & np.isfinite(df["value"])]


# ----------------------------------------------
# Fit AMiGA per (species, media, cond, replicate)
# ----------------------------------------------
records = []
group_cols = ["species", "media", "cond", "replicate"]

for key, sub in df.groupby(group_cols):
    species, media, cond, rep = key
    sub = sub.sort_values("time")

    t = sub["time"].values
    od = sub["value"].values

    try:
        m = GrowthModel(t, od)
        m.run()
        m.estimate_parameters()
        m.estimate_derivatives()

        # Extract primary values
        mumax = float(m.mu_max) if m.mu_max is not None else np.nan
        lag_time = float(m.lag_time) if m.lag_time is not None else np.nan

        # Derivative curve
        if m.derivative is None:
            max_der = np.nan
        else:
            max_der = float(np.nanmax(m.derivative[:, 1]))

        # Growth/no growth classification
        growth_call = int(m.has_growth())

    except Exception as e:
        print(f"AMiGA fitting failed for {key}: {e}", file=sys.stderr)
        mumax = np.nan
        lag_time = np.nan
        max_der = np.nan
        growth_call = 0

    records.append({
        "species": species,
        "media": media,
        "cond": cond,
        "replicate": rep,
        "mumax": mumax,
        "lag": lag_time,
        "max_derivative": max_der,
        "growth_call": growth_call
    })


out = pd.DataFrame(records)
os.makedirs(os.path.dirname(output_csv), exist_ok=True)
out.to_csv(output_csv, index=False)

print("Wrote:", output_csv)
