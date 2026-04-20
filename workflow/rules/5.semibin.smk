"""
Metagenome binning with Semibin
Metagenomic Binning Using Siamese Neural Networks for short and long reads

Here I am using this binning tool in single-sample binning mode, that is each assembly with its own bam file
Running CPU version
"""
from glob import glob
print(f"sample_names: {len(sample_names)}")

rule semibin_multi_sample_simka:
    input:
        assembly = os.path.join(dir_assembly,"{sample}.megahit.contigs.fa"),
        bam = os.path.join(dir_binning, "{sample}_cluster_50", "done.txt"),
    output:
        bins = os.path.join(dir_binning, "{sample}_semibin_bins", "done.txt")
    params:
        bin_dir= os.path.join(dir_binning, "{sample}_semibin_bins"),
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
        #generating the bins with semibin2
        SemiBin2 single_easy_bin -i {input.assembly} \
            -b {params.bam_dir}/*.bam -o {params.bin_dir} -t {threads}
        touch {output.bins}
        """
 
rule merge_semibins:
    input:
        ins=expand(bins = os.path.join(dir_binning, "{sample}_semibin_bins", "done.txt"), sample=sample_names)
    output:
        os.path.join(dir_binning, "merged_semibins", "done.txt")
    localrule: True
    params:
        dirs=os.path.join(dir_binning, "merged_semibins"),
        bins=os.path.join(dir_binning, "{sample}_semibin_bins", "output_bins"),
        samples="{sample}"
    shell:
        """
        #merging the bins from all samples into one directory
        mkdir -p {params.dirs}
        for sample in {params.samples}; do
            [ -e "$f" ] || continue
            bn=$(basename "$f")
            newname="${{sample}}_${{bn}}"
            cp -r {params.bins}/*.fa.gz {params.dirs}/$newname"
        done
        touch {output}
        """