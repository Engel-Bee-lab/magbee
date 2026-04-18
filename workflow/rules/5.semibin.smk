"""
Metagenome binning with Semibin
Metagenomic Binning Using Siamese Neural Networks for short and long reads

Here I am using this binning tool in single-sample binning mode, that is each assembly with its own bam file
Running CPU version
"""
rule semibin_multi_sample_simka:
    input:
        assembly = os.path.join(dir_assembly,"{sample}.megahit.contigs.fa"),
        bam = os.path.join(dir_binning, "{sample}_cluster_50", "done.txt"),
    output:
        bins = os.path.join(dir_binning, "{sample}_semibin_bins", "done.txt")
    params:
        bin_dir= os.path.join(dir_binning, "{sample}_semibin_bins"),
        temp= os.path.join(dir_binning, "{sample}_semibin_bins", "temp"),
        bam_dir= os.path.join(dir_binning, "{sample}_cluster_50"),
        sample="{sample}",
    conda:
        os.path.join(dir_env, "semibin2.yaml")
    resources:
        mem_mb =config['resources']['long_shortjob']['mem_mb'],
        runtime = config['resources']['long_shortjob']['runtime']
    threads:
        config['resources']['long_shortjob']['threads']
    shell:
        """
        #generating a concatenated fasta file for semibin2 but with only one assembly file 
        #This is being done for renaming the contigs 
        SemiBin2 concatenate_fasta -i {input.assembly} -o {params.temp}/{params.sample}_concatenated.fa

        #generating the bins with semibin2
        SemiBin2 multi_easy_bin -i {params.temp}/{params.sample}_concatenated.fa \
            -b {params.bam_dir}/*.bam -o {params.bin_dir} -t {threads}
        touch {output.bins}
        """
 