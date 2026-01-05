
"""
DRep dereplication of metabat2 bins 
"""
rule drep_metabat2_bins:
    input:
        bins_done = os.path.join(dir_binning, "all_metabat2_bins", "done.txt")
    output:
        drep_dir = os.path.join(dir_binning, "drep_metabat2_bins", "done.txt")
    params:
        outdir = os.path.join(dir_binning, "drep_metabat2_bins"),
        genomes_dir = os.path.join(dir_binning, "drep_metabat2_bins", "genomes"),
        bins_dir = os.path.join(dir_binning, "all_metabat2_bins")
    conda:
        os.path.join(dir_env, "drep.yaml")
    resources:
        mem  = config['resources']['bigjob']['mem'],
        time = config['resources']['bigjob']['time']
    threads:
        config['resources']['bigjob']['cpu']
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}
        # run dRep on the copied genomes
        dRep dereplicate {params.outdir} -g {params.bins_dir}/*.fa -comp 0 -con 1000 --clusterAlg average -p {threads}
        touch {output.drep_dir}
        """
