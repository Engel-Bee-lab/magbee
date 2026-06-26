"""
This rule is in progress. 
My idea is to add this in once I have a output I am happy with. 
"""

rule gtdbtk_dastool_individual:
    input:
        bins_dir = os.path.join(dir_binning, "das_tool", "dastool_DASTool_summary.tsv")
    output:
        gtdbtk_dir = os.path.join(dir_binning, "gtdbtk_output_dastool", "done.txt"),
        gtdbtk_summary = os.path.join(dir_binning, "gtdbtk_output_dastool", "classify", "gtdbtk.bac120.summary.tsv"),
        gtdbtk_ar_sumamry = os.path.join(dir_binning, "gtdbtk_output_dastool", "classify", "gtdbtk.ar53.summary.tsv")
    params:
        bins = os.path.join(dir_binning, "das_tool", "dastool_DASTool_bins"),
        outdir = os.path.join(dir_binning, "gtdbtk_dastool_output"),
        database = config["databases"]["gtdbtk_db"]
    conda:
        os.path.join(dir_env, "gtdbtk.yaml")
    resources:
        mem_mb  = config['resources']['assemblyjob']['mem_mb'],
        runtime = config['resources']['assemblyjob']['runtime']
    threads:
        config['resources']['assemblyjob']['threads']
    shell:
        """
        set -euo pipefail

        mkdir -p {params.outdir}
        export GTDBTK_DATA_PATH={params.database}
        gtdbtk identify --genome_dir {params.bins} --cpus {threads} --out_dir {params.outdir}/identify -x fa
        gtdbtk align --identify_dir {params.outdir}/identify --out_dir {params.outdir}/align --cpus {threads} 
        gtdbtk classify --genome_dir {params.bins} --out_dir {params.outdir}/classify --cpus {threads} -x fa -f --align_dir {params.outdir}/align
        touch {output.gtdbtk_dir}
        """

