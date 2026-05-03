"""
Metagenome binning with VAMB
Vamb is a family of metagenomic binners which feeds kmer composition and abundance into a variational autoencoder and clusters the embedding to form bins. 
Its binners perform excellently with multiple samples, and pretty good on single-sample data.
"""
from glob import glob

rule vamb_simka:
    input:
        contigs = os.path.join(dir_assembly,"{sample}.megahit.contigs.fa"),
        bam = os.path.join(dir_binning, "{sample}_cluster_50", "done.txt")
    output:
        abundance = os.path.join(dir_binning, "{sample}_vamb_bins", "vae_clusters_unsplit.tsv")
    params:
        contigs_rename = os.path.join(dir_binning, "{sample}_vamb_bins", "{sample}_contigs.fa"),
        bin_dir= os.path.join(dir_binning, "{sample}_vamb_bins"),
        bam_dir=os.path.join(dir_binning, "{sample}_cluster_50"),
    conda:
        os.path.join(dir_env, "vamb.yaml")
    resources:
        mem_mb =config['resources']['long_shortjob']['mem_mb'],
        runtime = config['resources']['long_shortjob']['runtime']
    threads:
        config['resources']['long_shortjob']['threads']
    shell:
        """
        #using the logic that one assembly was mapped with multiple samples from simka
        rm -rf {params.bin_dir}
        vamb bin default --outdir {params.bin_dir} --fasta {input.contigs} --bamdir {params.bam_dir} -t {threads}
        """

rule vamb_sep:
    input:
        contigs = os.path.join(dir_assembly,"{sample}.megahit.contigs.fa"),
        bam = os.path.join(dir_binning, "{sample}_vamb_bins", "vae_clusters_unsplit.tsv")
    output:
        bins = os.path.join(dir_binning, "{sample}_vamb_bins", "done.txt")
    params:
        bin_dir= os.path.join(dir_binning, "{sample}_vamb_bins"),
        outdir=os.path.join(dir_binning, "{sample}_vamb_bins", "bins"),
        scripts= os.path.join(dir_script, "vamb_bins_sep.py"),
    localrule: True
    shell:
        """
        python {params.scripts} --mapping {input.bam} --fasta {input.contigs} \
            --outdir {params.outdir}
        rm -rf {params.outdir}/bin_clustername.fasta
        touch {output.bins}
        """

rule vamb_bins:
    input:
        bins = expand(os.path.join(dir_binning, "{sample}_vamb_bins", "done.txt"), sample=sample_names)
    output:
        bins = os.path.join(dir_binning, "all_vamb_bins", "done.txt")
    localrule: True
    params:
        outdir= os.path.join(dir_binning, "all_vamb_bins"),
        sample=" ".join(sample_names)
    shell:
        """
        mkdir -p {params.outdir}
        for sample in {params.sample}; do
            src_dir={dir_binning}/${{sample}}_vamb_bins
            
            for f in "$src_dir"/bins/*.fna; do
                [ -e "$f" ] || continue
                bn=$(basename "$f")
                newname="${{sample}}_vamb_${{bn}}"
                cp "$f" "{params.outdir}/$newname"
            done
        done
        touch {output.bins}
        """