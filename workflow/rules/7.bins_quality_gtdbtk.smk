"""
GTDB on all the dastool bins
"""
rule gtdbtk_dastool_bins:
    input:
        bins_dir = os.path.join(dir_binning, "das_tool", "dastool_DASTool_summary.tsv")
    output:
        gtdbtk_summary = os.path.join(dir_binning, "gtdbtk_output_dastool", "classify", "gtdbtk.bac120.summary.tsv"),
        gtdbtk_ar_sumamry = os.path.join(dir_binning, "gtdbtk_output_dastool", "classify", "gtdbtk.ar53.summary.tsv")
    params:
        outdir = os.path.join(dir_binning, "gtdbtk_output_dastool"),
        database = config["databases"]["gtdbtk_db"],
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

        mkdir -p {params.outdir}
        export GTDBTK_DATA_PATH={params.database}

        gtdbtk identify --genome_dir "$BATCHDIR" --cpus {threads} --out_dir {params.outdir}/identify -x fa
        gtdbtk align --identify_dir {params.outdir}/identify --out_dir {params.outdir}/align --cpus {threads}
        gtdbtk classify --genome_dir "$BATCHDIR" --out_dir {params.outdir}/classify --cpus {threads} -x fa \
            --pplacer_cpus {threads} --align_dir {params.outdir}/align --debug

        touch {output.gtdbtk_summary}
        touch {output.gtdbtk_ar_sumamry}
        """

