"""
Rules to evaluate bins quality running Checkm2 
"""
from glob import glob

rule bins_checkm2:
    input:
        drep_dir = os.path.join(dir_binning, "drep_metabat2_bins", "done.txt")
    output:
        checkm2_dir = os.path.join(dir_binning, "drep_metabat2_bins", "checkm2_report.txt")
    params:
        bins=os.path.join(dir_binning, "{sample}_metabat2_bins"),
        outdir=os.path.join(dir_binning, "{sample}_metabat2_bins", "checkm2_output")
    conda:
        os.path.join(dir_env, "checkm2.yaml")
    resources:
        mem =config['resources']['smalljob']['mem'],
        time = config['resources']['smalljob']['time']
    threads: 
        config['resources']['smalljob']['cpu']
    shell:
        """
        mkdir -p {params.outdir}
        checkm2 assess -x fa -o 2 -t {threads} {params.bins} {params.outdir}
        cp {params.outdir}/checkm2_assessment.tsv {output.checkm2_dir}
        """