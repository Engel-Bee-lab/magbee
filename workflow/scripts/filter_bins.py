import argparse
import pandas as pd
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Filter bins into HQ and MQ based on completeness and contamination"
    )

    parser.add_argument("-i", "--input", required=True, help="Input TSV file (quality_report.tsv)")
    parser.add_argument("-o", "--outdir", default="filtered_bins", help="Output directory")

    parser.add_argument("--hq_comp", type=float, default=90.0, help="HQ completeness threshold (> )")
    parser.add_argument("--hq_cont", type=float, default=5.0, help="HQ contamination threshold (< )")

    parser.add_argument("--mq_comp_min", type=float, default=50.0, help="MQ min completeness (>= )")
    parser.add_argument("--mq_comp_max", type=float, default=90.0, help="MQ max completeness (<= )")
    parser.add_argument("--mq_cont", type=float, default=10.0, help="MQ contamination threshold (< )")

    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")

    # HQ bins
    hq = df[
        (df["Completeness"] > args.hq_comp) &
        (df["Contamination"] < args.hq_cont)
    ].copy()

    # MQ bins
    mq = df[
        (df["Completeness"] >= args.mq_comp_min) &
        (df["Completeness"] <= args.mq_comp_max) &
        (df["Contamination"] < args.mq_cont)
    ].copy()

    # remove HQ overlap from MQ
    mq = mq.drop(hq.index, errors="ignore")

    outdir = Path(args.outdir)
    hq_dir = outdir / "HQ"
    mq_dir = outdir / "MQ"

    hq_dir.mkdir(parents=True, exist_ok=True)
    mq_dir.mkdir(parents=True, exist_ok=True)

    # write outputs
    hq.to_csv(hq_dir / "HQ_bins.tsv", sep="\t", index=False)
    mq.to_csv(mq_dir / "MQ_bins.tsv", sep="\t", index=False)

    hq["Name"].to_csv(hq_dir / "HQ_bin_names.txt", index=False, header=False)
    mq["Name"].to_csv(mq_dir / "MQ_bin_names.txt", index=False, header=False)

    print(f"Done: {len(hq)} HQ bins, {len(mq)} MQ bins")
    print(f"Output written to: {outdir}")


if __name__ == "__main__":
    main()