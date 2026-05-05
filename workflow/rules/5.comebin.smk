"""
Binning rule using COMEBIN
"""
rule comebin_simka:
    input:
        assembly = os.path.join(dir_assembly,"{sample}.megahit.contigs.fa.gz"),
        bam=os.path.join(dir_backmapping, "{sample}_cluster_50", "done.txt")
    output:
        bins_dir = os.path.join(dir_binning, "{sample}_comebin_bins", "done.txt")
    params:
        outdir=os.path.join(dir_binning, "{sample}_comebin_bins"),
        bam_dir=os.path.join(dir_backmapping, "{sample}_cluster_50"),
        checkm_db = config['databases']['checkm_db']
    conda:
        os.path.join(dir_env, "comebin_gpu.yaml")
    resources:
        slurm_partition = "gpu",
        gres = "gpu:1",
        mem_mb =config['resources']['gpujob']['mem_mb'],
        runtime = config['resources']['gpujob']['runtime']
    threads:
        config['resources']['gpujob']['threads']
    shell:
        """
        rm -rf {params.outdir}
        checkm data setRoot {params.checkm_db}
        run_comebin.sh -a {input.assembly} -o {params.outdir} \
            -p {params.bam_dir} -t {threads}
        touch {output.bins_dir}
        """

rule comebin_all2all:
    input:
        assembly = os.path.join(dir_assembly,"{sample}.megahit.contigs.fa"),
        bam = os.path.join(dir_backmapping, "{sample}_bam", "done.txt")
    output:
        bins_dir = os.path.join(dir_binning, "{sample}_comebin_bins", "done.txt")
    params:
        outdir=os.path.join(dir_binning, "{sample}_comebin_bins"),
        bam_dir=os.path.join(dir_backmapping, "{sample}_bam"),
        checkm_db = config['databases']['checkm_db']
    conda:
        os.path.join(dir_env, "comebin_gpu.yaml")
    resources:
        slurm_partition = "gpu",
        gres = "gpu:1",
        mem_mb =config['resources']['gpujob']['mem_mb'],
        runtime = config['resources']['gpujob']['runtime']
    threads:
        config['resources']['gpujob']['threads']
    shell:
        """
        rm -rf {params.outdir}
        checkm data setRoot {params.checkm_db}
        run_comebin.sh -a {input.assembly} -o {params.outdir} \
            -p {params.bam_dir} -t {threads}
        touch {output.bins_dir}
        """

#this is work in progress still 
rule collect_comebin_bins:
    input:
        bins_done = expand(bins_dir = os.path.join(dir_binning, "{sample}_comebin_bins", "done.txt"), sample=sample_names),
    localrule: True
    output:
        collected_dir = os.path.join(dir_binning, "all_comebin_bins", "done.txt")
    params:
        outdir= os.path.join(dir_binning, "all_comebin_bins"),
        bin_dir= os.path.join(dir_binning, "{sample}_comebin_bins"),
        sample=" ".join(sample_names)
    shell:
        """
        mkdir -p {params.outdir}
        for sample in {params.sample}; do
            src_dir={params.bin_dir}/comebin_res/comebin_res_bins
            for f in "$src_dir"/*.fa*; do
                [ -e "$f" ] || continue
                bn=$(basename "$f")
                newname="${{sample}}_comebin_${{bn}}"
                cp "$f" "{params.outdir}/$newname"
            done
        done
        touch {output.collected_dir}
        """