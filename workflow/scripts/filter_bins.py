import argparse
import pandas as pd
from pathlib import Path
import shutil


def copy_bins(names, src_dir, dest_dir, ext=".fa"):
    dest_dir.mkdir(parents=True, exist_ok=True)

    for name in names:
        src_file = src_dir / f"{name}{ext}"

        if src_file.exists():
            shutil.copy2(src_file, dest_dir / src_file.name)
        else:
            print(f"[WARN] Missing file: {src_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Filter bins and copy FASTA files into ALL / HQ / MQ folders"
    )

    parser.add_argument("-i", "--input", required=True, help="quality_report.tsv")
    parser.add_argument("-o", "--outdir", default="filtered_bins")
    parser.add_argument("-s", "--src", required=True, help="Directory with bin FASTA files")
    parser.add_argument("--ext", default=".fa", help="FASTA extension (default .fa)")

    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")

    src_dir = Path(args.src)
    outdir = Path(args.outdir)

    # ALL / HQ / MQ logic
    all_bins = df.copy()

    hq = df[(df["Completeness"] > 90) & (df["Contamination"] < 5)].copy()

    mq = df[
        (df["Completeness"] >= 50) &
        (df["Completeness"] <= 90) &
        (df["Contamination"] < 10)
    ].copy()

    mq = mq.drop(hq.index, errors="ignore")

    # folders
    all_dir = outdir / "ALL"
    hq_dir = outdir / "HQ"
    mq_dir = outdir / "MQ"

    # copy files
    copy_bins(all_bins["Name"], src_dir, all_dir, args.ext)
    copy_bins(hq["Name"], src_dir, hq_dir, args.ext)
    copy_bins(mq["Name"], src_dir, mq_dir, args.ext)

    # also write TSV summaries
    all_bins.to_csv(all_dir / "all_bins.tsv", sep="\t", index=False)
    hq.to_csv(hq_dir / "HQ_bins.tsv", sep="\t", index=False)
    mq.to_csv(mq_dir / "MQ_bins.tsv", sep="\t", index=False)

    print(f"Done:")
    print(f"  ALL bins: {len(all_bins)}")
    print(f"  HQ bins : {len(hq)}")
    print(f"  MQ bins : {len(mq)}")
    print(f"Output: {outdir}")


if __name__ == "__main__":
    main()