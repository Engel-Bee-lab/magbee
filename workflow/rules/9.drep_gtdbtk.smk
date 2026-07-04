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

BATCH_SIZE = 50

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
        mem_mb = config['resources']['highmemjob']['mem_mb'],
        runtime = config['resources']['highmemjob']['runtime'],
    threads:
        config['resources']['highmemjob']['threads']
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
        gtdbtk classify --genome_dir "$BATCHDIR" --out_dir {params.outdir}/classify --cpus {threads} -x fasta \
            --pplacer_cpus {threads} --align_dir {params.outdir}/align --debug


        #rm -rf "$BATCHDIR"

        touch {output.bac_summary}
        touch {output.ar_summary}
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
        done = os.path.join(dir_species, "gtdbtk_output_derep", "done.txt"),
    run:
        # filter out zero-length files so the merge uses a real header source
        bac_files = [f for f in input.bac_summaries if os.path.getsize(f) > 0]
        ar_files = [f for f in input.ar_summaries if os.path.getsize(f) > 0]

        def merge(files, out_path):
            if not files:
                # create an empty file if no batch outputs exist
                os.makedirs(os.path.dirname(out_path), exist_ok=True)
                with open(out_path, "w") as out_f:
                    out_f.write("")
                return

            with open(out_path, "w") as out_f:
                with open(files[0]) as f:
                    out_f.write(f.read())
                for fname in files[1:]:
                    with open(fname) as f:
                        next(f, None)  # skip header
                        out_f.writelines(f)

        merge(bac_files, output.combined_bac)
        merge(ar_files, output.combined_ar)

        os.makedirs(os.path.dirname(output.done), exist_ok=True)
        with open(output.done, "w") as f:
            f.write("done\n")

