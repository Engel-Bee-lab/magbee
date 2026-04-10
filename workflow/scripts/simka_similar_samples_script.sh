#!/usr/bin/env bash
set -euo pipefail

# -----------------------------

# Parse arguments

# -----------------------------

usage() {
echo "Usage: $0 -s <sample> -d <dist.gz> -n <nr> -o <out_dir>"
exit 1
}

while getopts "s:d:n:o:" flag; do
case "${flag}" in
s) sample="${OPTARG}" ;;
d) dist="${OPTARG}" ;;
n) nr="${OPTARG}" ;;
o) out="${OPTARG}" ;;
*) usage ;;
esac
done

# -----------------------------

# Check inputs

# -----------------------------

[[ -z "${sample:-}" ]] && usage
[[ -z "${dist:-}" ]] && usage
[[ -z "${nr:-}" ]] && usage
[[ -z "${out:-}" ]] && usage

mkdir -p "$out"

outfile="${out}/${sample}_backmap_samples.txt"

echo "$sample" > "$outfile"

# -----------------------------

# Extract most similar samples

# -----------------------------

zcat "$dist" | awk -F';' -v sample="$sample" '
NR==1 {
for (i=2; i<=NF; i++) hdr[i]=$i
next
}
$1 == sample {
for (i=2; i<=NF; i++) {
if (hdr[i] != sample)
print hdr[i] ";" $i
}
}
' | sort -t';' -k2,2n 
| head -n "$nr" 
| cut -d';' -f1 
| sort >> "$outfile"


