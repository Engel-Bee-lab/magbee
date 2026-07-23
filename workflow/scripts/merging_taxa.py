#!/usr/bin/env python3
import argparse
import pandas as pd
import re


def _norm_taxon(x: str) -> str:
    if pd.isna(x):
        return ""
    x = str(x).strip()
    x = re.sub(r"\s+", " ", x)
    return x


def read_kraken_report(path: str) -> pd.DataFrame:
    # Kraken2 report columns:
    # 0:pct, 1:clade_reads, 2:taxon_reads, 3:rank_code, 4:taxid, 5:name
    k = pd.read_csv(path, sep="\t", header=None, names=[
        "kraken_pct", "kraken_clade_reads", "kraken_taxon_reads",
        "rank_code", "taxid", "taxon"
    ])
    k["taxon"] = k["taxon"].astype(str).str.strip()
    k["taxon_norm"] = k["taxon"].map(_norm_taxon)
    k["kraken_pct"] = pd.to_numeric(k["kraken_pct"], errors="coerce").fillna(0.0)
    k["kraken_frac"] = k["kraken_pct"] / 100.0
    return k[["taxid", "taxon", "taxon_norm", "rank_code", "kraken_pct", "kraken_frac", "kraken_clade_reads", "kraken_taxon_reads"]]


def read_rat_table(path: str) -> pd.DataFrame:
    # Flexible RAT parser: expects at least one taxonomy column + one abundance column
    r = pd.read_csv(path, sep="\t", comment="#")
    cols = {c.lower(): c for c in r.columns}

    taxon_col = None
    for key in ["lineage", "taxon", "name", "classification"]:
        if key in cols:
            taxon_col = cols[key]
            break
    if taxon_col is None:
        raise ValueError(f"No taxonomy column found in RAT file: {path}")

    abundance_col = None
    for key in ["abundance", "relative_abundance", "fraction", "percentage", "pct"]:
        if key in cols:
            abundance_col = cols[key]
            break
    if abundance_col is None:
        raise ValueError(f"No abundance column found in RAT file: {path}")

    out = pd.DataFrame({
        "taxon": r[taxon_col].astype(str).str.strip(),
        "rat_value_raw": pd.to_numeric(r[abundance_col], errors="coerce").fillna(0.0),
    })

    # If RAT values look like fraction (<=1), convert to pct too
    if out["rat_value_raw"].max() <= 1.0:
        out["rat_frac"] = out["rat_value_raw"]
        out["rat_pct"] = out["rat_value_raw"] * 100.0
    else:
        out["rat_pct"] = out["rat_value_raw"]
        out["rat_frac"] = out["rat_value_raw"] / 100.0

    out["taxon_norm"] = out["taxon"].map(_norm_taxon)
    return out[["taxon", "taxon_norm", "rat_pct", "rat_frac"]]


def merge_taxa(kraken_report: str, rat_table: str, sample: str) -> pd.DataFrame:
    k = read_kraken_report(kraken_report)
    r = read_rat_table(rat_table)

    m = pd.merge(k, r, on="taxon_norm", how="outer", suffixes=("_kraken", "_rat"))

    # Keep human-readable taxon
    m["taxon"] = m["taxon_kraken"].fillna(m["taxon_rat"]).fillna(m["taxon_norm"])

    m["kraken_pct"] = m["kraken_pct"].fillna(0.0)
    m["rat_pct"] = m["rat_pct"].fillna(0.0)
    m["kraken_frac"] = m["kraken_frac"].fillna(0.0)
    m["rat_frac"] = m["rat_frac"].fillna(0.0)

    # Simple consensus: mean of available non-zero estimates
    def consensus(row):
        vals = [v for v in [row["kraken_pct"], row["rat_pct"]] if v > 0]
        return sum(vals) / len(vals) if vals else 0.0

    m["consensus_pct"] = m.apply(consensus, axis=1)
    m["consensus_frac"] = m["consensus_pct"] / 100.0
    m["sample"] = sample

    return m[[
        "sample", "taxid", "rank_code", "taxon",
        "kraken_pct", "rat_pct", "consensus_pct",
        "kraken_frac", "rat_frac", "consensus_frac",
        "kraken_clade_reads", "kraken_taxon_reads"
    ]].sort_values("consensus_pct", ascending=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Merge Kraken2 + RAT abundance into one taxa table.")
    p.add_argument("--kraken-report", required=True, help="Kraken2 report file")
    p.add_argument("--rat-table", required=True, help="RAT abundance/classification table")
    p.add_argument("--sample", required=True, help="Sample name")
    p.add_argument("--out", required=True, help="Output TSV")
    args = p.parse_args()

    merged = merge_taxa(args.kraken_report, args.rat_table, args.sample)
    merged.to_csv(args.out, sep="\t", index=False)