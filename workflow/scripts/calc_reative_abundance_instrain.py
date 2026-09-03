#!/usr/bin/env python3
"""
Calculate relative abundance from inStrain genome_info.tsv output, without
qPCR normalization (Python port of the relevant chunk of Load_data.Rmd from
Aiswarya-prasad/honeybee-cross-species-metagenomics).

Two relative abundance metrics are computed, per sample:

1. Coverage-based relative abundance (rel_cov):
       rel_cov = coverage / sum(coverage across all genomes in that sample)

2. Read-count-based relative abundance (percentage_expected_number_of_genomes):
       expected_number_of_genomes = filtered_read_pair_count * read_length * 2 / genome_length
       percentage_expected_number_of_genomes = expected_number_of_genomes / sum(...) * 100
   (this approximates genome-copy-equivalents recruited per sample, then
   normalizes to a percentage -- an alternative to raw coverage that accounts
   for genome length differences)

Also reports the simpler read-share metric:
       percentage_reads = filtered_read_pair_count / sum(filtered_read_pair_count) * 100

Usage:
    python calc_relative_abundance.py All_genome_info.tsv rel_abundance_long.tsv \
        [--breadth_cutoff 0.5] [--read_length 150] [--matrix_out rel_abundance_matrix.tsv]

Input:
    All_genome_info.tsv must have a 'sample' column (as produced by the
    earlier awk command that prepends sample name) plus the standard inStrain
    genome_info.tsv columns: genome, coverage (or coverage_median), breadth,
    length, filtered_read_pair_count.

Output:
    rel_abundance_long.tsv - original table with rel_cov, percentage_reads,
        and percentage_expected_number_of_genomes columns added
    rel_abundance_matrix.tsv (optional, --matrix_out) - sample x genome
        matrix of rel_cov, NaN filled with 0 (breadth_cutoff applied)
"""
import argparse
import sys
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("genome_info_tsv")
    ap.add_argument("output_tsv")
    ap.add_argument("--breadth_cutoff", type=float, default=0.0,
                     help="minimum breadth to include a genome in relative abundance calcs (default: 0, i.e. no filter)")
    ap.add_argument("--read_length", type=float, default=150,
                     help="read length used for expected_number_of_genomes calc (default: 150)")
    ap.add_argument("--matrix_out", default=None,
                     help="optional path to write a sample x genome rel_cov matrix (wide format)")
    args = ap.parse_args()

    df = pd.read_csv(args.genome_info_tsv, sep="\t")

    required = {"sample", "genome", "breadth", "length", "filtered_read_pair_count"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"ERROR: missing required columns: {missing}. Available: {list(df.columns)}")

    # coverage column may be named 'coverage' or 'coverage_median' depending on inStrain version
    cov_col = "coverage" if "coverage" in df.columns else "coverage_median"
    if cov_col not in df.columns:
        sys.exit(f"ERROR: no coverage column found (looked for 'coverage' / 'coverage_median'). "
                  f"Available: {list(df.columns)}")

    # apply breadth filter (matches R: filter(breadth > breadth_cutoff))
    before = len(df)
    df = df[df["breadth"] > args.breadth_cutoff].copy()
    print(f"Kept {len(df)}/{before} rows after breadth > {args.breadth_cutoff} filter")

    # --- 1. coverage-based relative abundance ---
    df["coverage"] = df[cov_col]
    df["rel_cov"] = df.groupby("sample")["coverage"].transform(lambda x: x / x.sum())

    # --- 2. read-count-based relative abundance ---
    df["expected_number_of_genomes"] = (
        df["filtered_read_pair_count"] * args.read_length * 2 / df["length"]
    )
    df["percentage_expected_number_of_genomes"] = df.groupby("sample")["expected_number_of_genomes"] \
        .transform(lambda x: x / x.sum() * 100)
    df["percentage_reads"] = df.groupby("sample")["filtered_read_pair_count"] \
        .transform(lambda x: x / x.sum() * 100)

    df.to_csv(args.output_tsv, sep="\t", index=False)
    print(f"Wrote {len(df)} rows to {args.output_tsv}")

    if args.matrix_out:
        matrix = (
            df[["sample", "genome", "rel_cov"]]
            .drop_duplicates()
            .pivot(index="sample", columns="genome", values="rel_cov")
            .fillna(0)
        )
        matrix.to_csv(args.matrix_out, sep="\t")
        print(f"Wrote sample x genome rel_cov matrix ({matrix.shape[0]}x{matrix.shape[1]}) to {args.matrix_out}")


if __name__ == "__main__":
    main()