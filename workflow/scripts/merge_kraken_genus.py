#!/usr/bin/env python3
"""
Merge per-sample Kraken2 report files into one genus-level wide table.

Standard Kraken2 report columns (tab-separated, no header):
    1. percent          - % of reads covered by the clade rooted at this taxon
    2. clade_reads      - # reads in the clade rooted at this taxon (this taxon + everything below it)
    3. taxon_reads      - # reads assigned directly to this taxon (not resolved further)
    4. rank_code        - U, R, D, K, P, C, O, F, G, S, (sometimes G1/G2/... for sub-genus groupings)
    5. taxid            - NCBI taxonomy ID
    6. name             - indented scientific name

For genus-level abundance you want `clade_reads` on rows where rank_code == "G": that
value already includes every species (and strain) classified beneath that genus, so it
is the correct "total reads belonging to this genus" number -- no further summing needed.

Output table:
    rows    = genus name, kept together with its taxid
    columns = for each metagenome <sample>:
                  <sample>_reads    (clade_reads at the genus-level row)
                  <sample>_percent  (percent at the genus-level row)

Usage:
    python merge_kraken_genus.py --input-dir /path/to/kraken_reports --pattern "*.kreport" \
        --output merged_kraken_genus.tsv

    By default only rows with rank_code exactly "G" are kept. If your database reports
    sub-genus ranks (G1, G2, ...) and you want those folded in too, pass --rank-prefix G
    instead of the default exact match (see --rank / --rank-prefix below).
"""

import argparse
import glob
import os
import sys

import pandas as pd

KRAKEN_REPORT_COLS = ["percent", "clade_reads", "taxon_reads", "rank_code", "taxid", "name"]


def sample_name_from_path(path: str, strip_suffixes) -> str:
    base = os.path.basename(path)
    for suf in strip_suffixes:
        if base.endswith(suf):
            base = base[: -len(suf)]
            break
    else:
        base = os.path.splitext(base)[0]
    return base


def load_one_kraken_report(path: str, sample: str, rank: str, rank_prefix: bool) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=KRAKEN_REPORT_COLS,
        dtype={"percent": float, "clade_reads": "Int64", "taxon_reads": "Int64", "rank_code": str, "taxid": "Int64"},
    )
    if rank_prefix:
        df = df[df["rank_code"].str.startswith(rank, na=False)].copy()
    else:
        df = df[df["rank_code"] == rank].copy()

    df["name"] = df["name"].str.strip()
    df = df[["name", "taxid", "clade_reads", "percent"]].copy()
    df = df.rename(
        columns={
            "clade_reads": f"{sample}_reads",
            "percent": f"{sample}_percent",
        }
    )
    return df


def main():
    ap = argparse.ArgumentParser(description="Merge per-sample Kraken2 reports into one genus-level wide table.")
    input_group = ap.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--input-dir",
        help="Directory containing one Kraken2 report file per metagenome (used with --pattern).",
    )
    input_group.add_argument(
        "--inputs",
        nargs="+",
        help=(
            "Explicit list of Kraken2 report file paths (space-separated), e.g. the resolved "
            "output of a Snakemake expand() over per-sample subdirectories. Use this instead of "
            "--input-dir when your reports aren't all in one flat directory."
        ),
    )
    ap.add_argument(
        "--pattern",
        default="*.kreport",
        help="Glob pattern (relative to --input-dir) matching Kraken2 report files. Only used with --input-dir. Default: *.kreport",
    )
    ap.add_argument(
        "--rank",
        default="G",
        help="Rank code to keep. Default: 'G' (genus). Use 'S' for species, 'F' for family, etc.",
    )
    ap.add_argument(
        "--rank-prefix",
        action="store_true",
        help="Match rank codes by prefix (e.g. 'G' also matches 'G1', 'G2' sub-genus rows) instead of exact match.",
    )
    ap.add_argument(
        "--strip-suffix",
        action="append",
        default=None,
        help=(
            "Suffix to strip from each filename to derive the sample/metagenome name "
            "(can be given multiple times, first match wins). "
            "Default tries '.kreport', '_kraken2_report.txt', '.report.txt', '.txt' in that order."
        ),
    )
    ap.add_argument("--output", required=True, help="Path for the merged output table (.tsv or .csv).")
    ap.add_argument("--sep", default="\t", help="Field separator for the output file. Default: tab.")
    args = ap.parse_args()

    strip_suffixes = args.strip_suffix or [
        ".kreport",
        ".kraken2.report.txt",
        "_kraken2_report.txt",
        ".report.txt",
        ".txt",
    ]

    if args.inputs:
        files = sorted(args.inputs)
        missing = [f for f in files if not os.path.isfile(f)]
        if missing:
            sys.exit(f"These --inputs paths don't exist: {missing}")
    else:
        files = sorted(glob.glob(os.path.join(args.input_dir, args.pattern)))
        if not files:
            sys.exit(f"No files matched {os.path.join(args.input_dir, args.pattern)!r}")

    print(f"Found {len(files)} Kraken2 report files.")

    merged = None
    seen_samples = set()
    for path in files:
        sample = sample_name_from_path(path, strip_suffixes)
        if sample in seen_samples:
            sys.exit(f"Duplicate sample name '{sample}' derived from {path} — adjust --strip-suffix.")
        seen_samples.add(sample)

        df = load_one_kraken_report(path, sample, args.rank, args.rank_prefix)

        if merged is None:
            merged = df
        else:
            merged = merged.merge(df, on=["name", "taxid"], how="outer")

    read_cols = [c for c in merged.columns if c.endswith("_reads")]
    pct_cols = [c for c in merged.columns if c.endswith("_percent")]
    merged[read_cols] = merged[read_cols].fillna(0)
    merged[pct_cols] = merged[pct_cols].fillna(0.0)
    merged["taxid"] = merged["taxid"].astype("Int64")

    ordered_cols = ["name", "taxid"]
    for path in files:
        sample = sample_name_from_path(path, strip_suffixes)
        ordered_cols += [f"{sample}_reads", f"{sample}_percent"]
    merged = merged[ordered_cols]

    merged = merged.set_index("name")
    merged = merged.sort_index()

    merged.to_csv(args.output, sep=args.sep)
    print(f"Wrote merged table: {args.output}  ({merged.shape[0]} genera x {len(files)} metagenomes)")


if __name__ == "__main__":
    main()