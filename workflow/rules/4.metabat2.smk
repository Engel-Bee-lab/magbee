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
        mem =config['resources']['bigjob']['mem'],
        time = config['resources']['bigjob']['time']
    threads: 
        config['resources']['bigjob']['cpu']
    shell:
        """
        mkdir -p {params.outdir}
        jgi_summarize_bam_contig_depths --outputDepth {params.outdir}/depth.txt {params.bam_dir}/*.bam
        metabat2 -i {input.assembly} -a {params.outdir}/depth.txt -o {params.outdir}/bin -t {threads}
        touch {output.bins_dir}
        """