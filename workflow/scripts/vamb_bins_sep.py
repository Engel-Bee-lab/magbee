import argparse
from collections import defaultdict
from Bio import SeqIO
import os
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("--mapping", required=True)
parser.add_argument("--fasta", required=True)
parser.add_argument("--outdir", required=True)
parser.add_argument(
    "--min_size",
    type=int,
    default=200000,
    help="Minimum total bin length (bp) required to write a bin [default: 200000]"
)
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

clusters = defaultdict(list)
with open(args.mapping) as f:
    for line in f:
        cluster, contig = line.strip().split()
        clusters[cluster].append(contig)

if args.fasta.endswith(".gz"):
    handle = gzip.open(args.fasta, "rt")
else:
    handle = open(args.fasta, "r")

fasta_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
handle.close()

for cluster, contigs in clusters.items():

    total_len = sum(
        len(fasta_dict[c].seq)
        for c in contigs
        if c in fasta_dict
    )

    if total_len < args.min_size:
        print(
            f"Skipping cluster {cluster}: "
            f"{total_len:,} bp < {args.min_size:,} bp"
        )
        continue

    outfile = os.path.join(args.outdir, f"bin_{cluster}.fasta")
    with open(outfile, "w") as out:
        for c in contigs:
            if c in fasta_dict:
                SeqIO.write(fasta_dict[c], out, "fasta")

    print(
        f"Wrote {outfile} "
        f"({len(contigs)} contigs, {total_len:,} bp)"
    )