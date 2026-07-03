"""
Moved gtdbtk here after dereplication to avoid running gtdbtk on redundant bins. Also pplacer 
doesnt scale well with many bins, so this is a good way to reduce the number of bins that need to be placed.

This rule does run bins in batches if there are more than 100 bins even after derep 
"""
import os
import glob
import math

import os
import glob
import math

BATCH_SIZE = 100

derep_bins_dir = os.path.join(dir_species, "drep_dastools", "dereplicated_genomes")

# Computed once, at parse time, from whatever is on disk right now.
# If derep_bins_dir doesn't exist yet, this falls back to 1 batch so the
# DAG can still build -- you'll need to re-run snakemake after das_tool
# finishes to pick up the real file count.
if os.path.isdir(derep_bins_dir):
    _fa_files = sorted(glob.glob(os.path.join(derep_bins_dir, "*.fa*")))
    NUM_BATCHES = max(1, math.ceil(len(_fa_files) / BATCH_SIZE))
else:
    NUM_BATCHES = 1

def get_batch_files(batch_num, batch_size=BATCH_SIZE):
    fa_files = sorted(glob.glob(os.path.join(derep_bins_dir, "*.fa*")))
    start = batch_num * batch_size
    return fa_files[start:start + batch_size]


rule format_extn:
    input:
        summary = os.path.join(dir_species, "drep_dastools", "done.txt"),
    output:
        formatted = os.path.join(dir_temp, "reformatted_dastool_bin_extn.txt")
    params:
        derep_bins_dir = os.path.join(dir_species, "drep_dastools", "dereplicated_genomes")
    localrule: True
    shell:
        """
        set -euo pipefail

        for f in {params.derep_bins_dir}/*; do
            case "$f" in
                *.fasta.gz)
                    newname="${{f%.fasta.gz}}.fa"
                    gunzip -c "$f" > "$newname"
                    rm "$f"
                    ;;
                *.fa.gz)
                    newname="${{f%.fa.gz}}.fa"
                    gunzip -c "$f" > "$newname"
                    rm "$f"
                    ;;
                *.fasta)
                    mv "$f" "${{f%.fasta}}.fa"
                    ;;
                *.fa)
                    : # already correct
                    ;;
            esac
        done

        touch {output.formatted}
        """

rule gtdbtk_dastool_batch:
    input:
        summary = os.path.join(dir_species, "drep_dastools", "done.txt"),
    output:
        bac_summary = os.path.join(dir_species, "gtdbtk_output_derep", "batch_{batch_num}", "classify", "gtdbtk.bac120.summary.tsv"),
        ar_summary = os.path.join(dir_species, "gtdbtk_output_derep", "batch_{batch_num}", "classify", "gtdbtk.ar53.summary.tsv"),
    params:
        outdir = lambda wc: os.path.join(dir_species, "gtdbtk_output_derep", f"batch_{wc.batch_num}"),
        database = config["databases"]["gtdbtk_db"],
        batch_files = lambda wc: get_batch_files(int(wc.batch_num)),
    conda:
        os.path.join(dir_env, "gtdbtk.yaml")
    resources:
        mem_mb = config['resources']['assemblyjob']['mem_mb'],
        runtime = config['resources']['assemblyjob']['runtime'],
    threads:
        config['resources']['assemblyjob']['threads']
    shell:
        """
        set -euo pipefail

        BATCHDIR=batch_temp_{wildcards.batch_num}
        rm -rf "$BATCHDIR"
        mkdir -p "$BATCHDIR"

        for f in {params.batch_files}; do
            ln -s "$(readlink -f "$f")" "$BATCHDIR/$(basename "$f")"
        done

        mkdir -p {params.outdir}
        export GTDBTK_DATA_PATH={params.database}

        gtdbtk identify --genome_dir "$BATCHDIR" --cpus {threads} --out_dir {params.outdir}/identify -x fasta
        gtdbtk align --identify_dir {params.outdir}/identify --out_dir {params.outdir}/align --cpus {threads}
        gtdbtk classify --genome_dir "$BATCHDIR" --out_dir {params.outdir}/classify --cpus {threads} -x fasta -f --align_dir {params.outdir}/align

        rm -rf "$BATCHDIR"
        """


rule combine_gtdbtk_results:
    input:
        bac_summaries = expand(
            os.path.join(dir_species, "gtdbtk_output_derep", "batch_{batch_num}", "classify", "gtdbtk.bac120.summary.tsv"),
            batch_num=range(NUM_BATCHES),
        ),
        ar_summaries = expand(
            os.path.join(dir_species, "gtdbtk_output_derep", "batch_{batch_num}", "classify", "gtdbtk.ar53.summary.tsv"),
            batch_num=range(NUM_BATCHES),
        ),
    output:
        combined_bac = os.path.join(dir_reports, "gtdbtk_output_derep", "gtdbtk.bac120.summary.tsv"),
        combined_ar = os.path.join(dir_reports, "gtdbtk_output_derep", "gtdbtk.ar53.summary.tsv"),
    run:
        def merge(files, out_path):
            with open(out_path, "w") as out_f:
                with open(files[0]) as f:
                    out_f.write(f.read())
                for fname in files[1:]:
                    with open(fname) as f:
                        next(f, None)  # skip header
                        out_f.writelines(f)

        merge(list(input.bac_summaries), output.combined_bac)
        merge(list(input.ar_summaries), output.combined_ar)

        with open(output.done, "w") as f:
            f.write("done\n")

