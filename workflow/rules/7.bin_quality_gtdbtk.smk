"""
This rule is in progress. 
My idea is to add this in once I have a output I am happy with. 
"""

rule gtdbtk_bins:
    input:
        bins_done = os.path.join(dir_binning, "all_metabat2_bins", "done.txt"),
    output:
        gtdbtk_dir = os.path.join(dir_binning, "gtdbtk_output", "done.txt")
    params:
        bins = os.path.join(dir_binning, "all_metabat2_bins"),
        outdir = os.path.join(dir_binning, "gtdbtk_output"),
        database = config["databases"]["gtdbtk_db"]
    conda:
        os.path.join(dir_env, "gtdbtk.yaml")
    resources:
        mem_mb  = config['resources']['bigjob']['mem_mb'],
        runtime = config['resources']['bigjob']['runtime']
    threads:
        config['resources']['bigjob']['threads']
    shell:
        """
        set -euo pipefail

        mkdir -p {params.outdir}
        export GTDBTK_DATA_PATH={params.database}
        gtdbtk identify --genome_dir {params.bins} --cpus {threads} --out_dir {params.outdir} -x fa
        gtdbtk align --identify_dir {params.outdir} --out_dir {params.outdir} --cpus {threads} 
        gtdbtk classify --genome_dir {params.bins} --out_dir {params.outdir} --cpus {threads} -x fa -f --align_dir {params.outdir}
        touch {output.gtdbtk_dir}
        """

