import argparse
from collections import defaultdict
from Bio import SeqIO
import os
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("--mapping", required=True)
parser.add_argument("--fasta", required=True)
parser.add_argument("--outdir", required=True)
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
    
fasta_dict = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))

for cluster, contigs in clusters.items():
    outfile = os.path.join(args.outdir, f"bin_{cluster}.fasta")
    with open(outfile, "w") as out:
        for c in contigs:
            if c in fasta_dict:
                SeqIO.write(fasta_dict[c], out, "fasta")