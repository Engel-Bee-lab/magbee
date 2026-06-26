
"""
DRep dereplication of metabat2 bins 
"""
rule drep_dastool_bins:
    input:
        quality=os.path.join(dir_reports, "quality_mags.done")
    output:
        drep_dir = os.path.join(dir_species, "drep_dastool", "done.txt")
        samples= os.path.join(dir_species, "{mag}_bams", "{sample}.bam"),
    params:
        quality_mags=os.path.join(dir_reports, "quality_mags"),
        genomes_dir = os.path.join(dir_species, "drep_dastools"),
    conda:
        os.path.join(dir_env, "drep.yaml")
    resources:
        mem_mb  = config['resources']['bigjob']['mem_mb'],
        runtime = config['resources']['bigjob']['runtime']
    threads:
        config['resources']['bigjob']['threads']
    shell:
        """
        set -euo pipefail
        mkdir -p {params.genomes_dir}
        # run dRep on the copied genomes
        dRep dereplicate {params.genomes_dir} -g {params.quality_mags}/*.fa -comp 0 -con 1000 --clusterAlg average -p {threads}
        touch {output.drep_dir}
        """
