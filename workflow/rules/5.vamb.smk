"""
Metagenome binning with VAMB
Vamb is a family of metagenomic binners which feeds kmer composition and abundance into a variational autoencoder and clusters the embedding to form bins. 
Its binners perform excellently with multiple samples, and pretty good on single-sample data.
"""

rule vamb_simka:
    input:
        contigs = os.path.join(dir_assembly,"{sample}.megahit.contigs.fa"),
        bam=os.path.join(dir_binning, "{sample}_cluster_50", "done.txt")
    output:
        bins = os.path.join(dir_binning, "{sample}_vamb_bins", "done.txt")
    params:
        bin_dir= os.path.join(dir_binning, "{sample}_vamb_bins"),
        bam_dir=os.path.join(dir_binning, "{sample}_cluster_50")
    conda:
        os.path.join(dir_env, "vamb.yaml")
    resources:
        mem_mb =config['resources']['bigjob']['mem_mb'],
        runtime = config['resources']['bigjob']['runtime']
    threads:
        config['resources']['bigjob']['threads']
    shell:
        """
        #using the logic that one assembly was mapped with multiple samples from simka

        vamb -o C -p {threads} --outdir {params.bin_dir} --fasta {input.contigs} --bamfiles {params.bam_dir}/*.bam
        touch {output.bins}
        """
