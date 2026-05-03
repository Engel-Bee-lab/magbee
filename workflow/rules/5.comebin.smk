rule comebin_simka:
    input:
        assembly = os.path.join(dir_assembly,"{sample}.megahit.contigs.fa"),
        bam=os.path.join(dir_binning, "{sample}_cluster_50", "done.txt")
    output:
        bins_dir = os.path.join(dir_binning, "{sample}_comebin_bins", "done.txt")
    params:
        outdir=os.path.join(dir_binning, "{sample}_comebin_bins"),
        bam_dir=os.path.join(dir_binning, "{sample}_cluster_50")
        checkm_db = config['databases']['checkm']
    conda:
        os.path.join(dir_env, "comebin.yaml")
    resources:
        mem_mb =config['resources']['bigjob']['mem_mb'],
        runtime = config['resources']['bigjob']['runtime']
    threads:
        config['resources']['bigjob']['threads']
    shell:
        """
        mkdir -p {params.outdir}
        checkm data setRoot {params.checkm_db}
        run_comebin.sh -a {input.assembly} -o {params.outdir} \
            -p {params.bam_dir}/*.bam -t {threads}
        touch {output.bins_dir}
        """

#this is work in progress still 
"""
rule collect_comebin_bins:
    input:
        bins_done = expand(bins_dir = os.path.join(dir_binning, "{sample}_comebin_bins", "done.txt"), sample=sample_names),
    localrule: True
    output:
        collected_dir = os.path.join(dir_binning, "all_comebin_bins", "done.txt")
    params:
        outdir= os.path.join(dir_binning, "all_comebin_bins"),
        sample=" ".join(sample_names)
    shell:
        mkdir -p {params.outdir}
        for sample in {params.sample}; do
            src_dir={dir_binning}/${{sample}}_metabat2_bins
            for f in "$src_dir"/bin*.fa; do
                [ -e "$f" ] || continue
                bn=$(basename "$f")
                newname="${{sample}}_${{bn}}"
                cp "$f" "{params.outdir}/$newname"
            done
        done
        touch {output.collected_dir}
"""