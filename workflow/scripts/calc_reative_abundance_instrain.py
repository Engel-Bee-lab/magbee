#!/usr/bin/env python3
"""
Calculate relative abundance from inStrain genome_info.tsv output, without
qPCR normalization (Python port of the relevant chunk of Load_data.Rmd from
Aiswarya-prasad/honeybee-cross-species-metagenomics).

IMPORTANT (matches the paper's actual order of operations): rel_cov and the
other per-sample metrics are computed on the UNFILTERED data first -- the
denominator is the sum across every genome inStrain reported for that
sample, regardless of breadth. The breadth > breadth_cutoff filter is only
applied afterwards, as a "detected" flag, and when building the final
sample x genome matrix. Filtering before computing rel_cov (as an earlier
version of this script did) changes the denominator and does not match the
paper's code.

Two relative abundance metrics are computed, per sample:

1. Coverage-based relative abundance (rel_cov):
       rel_cov = coverage / sum(coverage across ALL genomes in that sample)

2. Read-count-based relative abundance (percentage_expected_number_of_genomes):
       expected_number_of_genomes = filtered_read_pair_count * read_length * 2 / genome_length
       percentage_expected_number_of_genomes = expected_number_of_genomes / sum(...) * 100
   (this approximates genome-copy-equivalents recruited per sample, then
   normalizes to a percentage -- an alternative to raw coverage that accounts
   for genome length differences)

Also reports the simpler read-share metric:
       percentage_reads = filtered_read_pair_count / sum(filtered_read_pair_count) * 100

And a boolean "detected" flag (breadth > breadth_cutoff), matching the
paper's detected-genome logic -- used to build the final matrix, but NOT
used when computing the rel_cov/percentage columns above.

Usage:
    # single pre-merged table (already has a 'sample' column):
    python calc_relative_abundance.py All_genome_info.tsv rel_abundance_long.tsv \
        [--breadth_cutoff 0.5] [--read_length 150] [--matrix_out rel_abundance_matrix.tsv]
        [--prevalence_out prevalence.tsv]

    # OR merge a batch of inStrain profile directories directly (one call does both steps):
    python calc_relative_abundance.py All_genome_info.tsv rel_abundance_long.tsv \
        --profile_glob "instrain_profiles/*/output/*_genome_info.tsv" \
        --matrix_out rel_abundance_matrix.tsv --prevalence_out prevalence.tsv

Input:
    Either:
    (a) All_genome_info.tsv already merged, with a 'sample' column (as
        produced by the earlier awk command that prepends sample name), or
    (b) --profile_glob pointing at each sample's own genome_info.tsv under
        its inStrain profile output directory (e.g. what
        "inStrain profile ... -o <sample>" produces at
        <sample>/output/<sample>_genome_info.tsv). The sample name is taken
        from the filename itself (stripping the _genome_info.tsv suffix),
        so it works regardless of directory nesting. The merged table is
        written out to the genome_info_tsv path given as the first
        argument, so you get it as an artifact too.
    Either way, the table needs the standard inStrain genome_info.tsv
    columns: genome, coverage (or coverage_median), breadth, length,
    filtered_read_pair_count.

Output:
    rel_abundance_long.tsv - full (unfiltered) table with rel_cov,
        percentage_reads, percentage_expected_number_of_genomes, and
        detected columns added
    rel_abundance_matrix.tsv (optional, --matrix_out) - sample x genome
        matrix of rel_cov, restricted to detected rows (breadth > cutoff),
        NaN filled with 0
    prevalence.tsv (optional, --prevalence_out) - per-genome prevalence,
        i.e. fraction of samples where it was detected
"""
import argparse
import glob
import os
import sys
import pandas as pd


def load_and_merge_profiles(profile_glob):
    """Read every sample's genome_info.tsv matched by profile_glob, tag each
    with its sample name (from the filename), and concatenate."""
    files = sorted(glob.glob(profile_glob))
    if not files:
        sys.exit(f"ERROR: --profile_glob '{profile_glob}' matched no files")
    dfs = []
    for f in files:
        sample = os.path.basename(f)
        for suffix in ("_genome_info.tsv", ".tsv"):
            if sample.endswith(suffix):
                sample = sample[: -len(suffix)]
                break
        d = pd.read_csv(f, sep="\t")
        d.insert(0, "sample", sample)
        dfs.append(d)
    print(f"Merged {len(files)} profile genome_info.tsv files ({sum(len(d) for d in dfs)} rows total)")
    return pd.concat(dfs, ignore_index=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("genome_info_tsv",
                     help="path to a pre-merged table (input), OR, with --profile_glob, "
                          "the path the merged table will be WRITTEN to")
    ap.add_argument("output_tsv")
    ap.add_argument("--profile_glob", default=None,
                     help="glob pattern matching each sample's genome_info.tsv under its inStrain "
                          "profile output dir, e.g. 'instrain_profiles/*/output/*_genome_info.tsv' "
                          "-- if given, genome_info_tsv is treated as an output path for the merged table")
    ap.add_argument("--breadth_cutoff", type=float, default=0.5,
                     help="minimum breadth for a genome to be 'detected' (default: 0.5, matches the paper)")
    ap.add_argument("--read_length", type=float, default=150,
                     help="read length used for expected_number_of_genomes calc (default: 150)")
    ap.add_argument("--matrix_out", default=None,
                     help="optional path to write a sample x genome rel_cov matrix (wide format, detected-only)")
    ap.add_argument("--prevalence_out", default=None,
                     help="optional path to write per-genome prevalence (fraction of samples detected in)")
    args = ap.parse_args()

    if args.profile_glob:
        df = load_and_merge_profiles(args.profile_glob)
        df.to_csv(args.genome_info_tsv, sep="\t", index=False)
        print(f"Wrote merged table ({len(df)} rows) to {args.genome_info_tsv}")
    else:
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

    # --- compute on the FULL, unfiltered table first (matches Load_data.Rmd order) ---
    df["coverage"] = df[cov_col]
    df["rel_cov"] = df.groupby("sample")["coverage"].transform(lambda x: x / x.sum())

    df["expected_number_of_genomes"] = (
        df["filtered_read_pair_count"] * args.read_length * 2 / df["length"]
    )
    df["percentage_expected_number_of_genomes"] = df.groupby("sample")["expected_number_of_genomes"] \
        .transform(lambda x: x / x.sum() * 100)
    df["percentage_reads"] = df.groupby("sample")["filtered_read_pair_count"] \
        .transform(lambda x: x / x.sum() * 100)

    # "detected" flag -- applied AFTER the above, never before
    df["detected"] = df["breadth"] > args.breadth_cutoff
    print(f"{df['detected'].sum()}/{len(df)} rows flagged detected (breadth > {args.breadth_cutoff})")

    df.to_csv(args.output_tsv, sep="\t", index=False)
    print(f"Wrote {len(df)} rows to {args.output_tsv}")

    if args.matrix_out:
        matrix = (
            df[df["detected"]][["sample", "genome", "rel_cov"]]
            .drop_duplicates()
            .pivot(index="sample", columns="genome", values="rel_cov")
            .fillna(0)
        )
        matrix.to_csv(args.matrix_out, sep="\t")
        print(f"Wrote sample x genome rel_cov matrix ({matrix.shape[0]}x{matrix.shape[1]}) to {args.matrix_out}")

    if args.prevalence_out:
        n_samples = df["sample"].nunique()
        prevalence = (
            df[df["detected"]].groupby("genome")["sample"].nunique() / n_samples
        ).sort_values(ascending=False)
        prevalence.to_csv(args.prevalence_out, sep="\t", header=["prevalence"])
        print(f"Wrote per-genome prevalence ({len(prevalence)} genomes, {n_samples} samples) to {args.prevalence_out}")


if __name__ == "__main__":
    main()