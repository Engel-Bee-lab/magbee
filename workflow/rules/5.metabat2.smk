"""
Binning rule using metabat2
"""
from glob import glob

rule metabat2_binning_individual_sample:
    input:
        assembly = os.path.join(dir_assembly,"{sample}.megahit.contigs.fa"),
        bam = os.path.join(dir_binning, "{sample}_bam", "done.txt"),
    output:
        bins_dir = os.path.join(dir_binning, "{sample}_metabat2_bins", "done.txt")
    params:
        outdir=os.path.join(dir_binning, "{sample}_metabat2_bins"),
        bam_dir=os.path.join(dir_binning, "{sample}_bam")
    conda:
        os.path.join(dir_env, "metabat2.yaml")
    resources:
        mem_mb =config['resources']['bigjob']['mem_mb'],
        runtime = config['resources']['bigjob']['runtime']
    threads: 
        config['resources']['bigjob']['threads']
    shell:
        """
        mkdir -p {params.outdir}
        jgi_summarize_bam_contig_depths --outputDepth {params.outdir}/depth.txt {params.bam_dir}/*.bam
        metabat2 -i {input.assembly} -a {params.outdir}/depth.txt -o {params.outdir}/bin -t {threads}
        touch {output.bins_dir}
        """

rule metabat2_binning_simka_samples:
    input:
        assembly = os.path.join(dir_assembly,"{sample}.megahit.contigs.fa"),
        bam=os.path.join(dir_binning, "{sample}_cluster_50", "done.txt")
    output:
        bins_dir = os.path.join(dir_binning, "{sample}_metabat2_bins", "done.txt")
    params:
        outdir=os.path.join(dir_binning, "{sample}_metabat2_bins"),
        bam_dir=os.path.join(dir_binning, "{sample}_cluster_50")
    conda:
        os.path.join(dir_env, "metabat2.yaml")
    resources:
        mem_mb =config['resources']['bigjob']['mem_mb'],
        runtime = config['resources']['bigjob']['runtime']
    threads:
        config['resources']['bigjob']['threads']
    shell:
        """
        mkdir -p {params.outdir}
        jgi_summarize_bam_contig_depths --outputDepth {params.outdir}/depth.txt {params.bam_dir}/*.bam
        metabat2 -i {input.assembly} -a {params.outdir}/depth.txt -o {params.outdir}/bin -t {threads}
        touch {output.bins_dir}
        """

rule collect_metabat2_bins:
    input:
        bins_done = expand(os.path.join(dir_binning, "{sample}_metabat2_bins", "done.txt"), sample=sample_names),
    localrule: True
    output:
        collected_dir = os.path.join(dir_binning, "all_metabat2_bins", "done.txt")
    params:
        outdir= os.path.join(dir_binning, "all_metabat2_bins"),
        sample=" ".join(sample_names)
    shell:
        """
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
