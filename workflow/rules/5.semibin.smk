"""
Metagenome binning with Semibin
Metagenomic Binning Using Siamese Neural Networks for short and long reads

Here I am using this binning tool in single-sample binning mode, that is each assembly with its own bam file
Running CPU version
"""
from glob import glob
#print(f"sample_names: {len(sample_names)})")

rule semibin_multi_sample_simka:
    input:
        assembly = os.path.join(dir_assembly,"{sample}.megahit.contigs.fa.gz"),
        bam = os.path.join(dir_backmapping, "{sample}_cluster_50", "done.txt")
    output:
        bins = os.path.join(dir_binning, "{sample}_semibin_bins", "done.txt")
    params:
        bin_dir= os.path.join(dir_binning, "{sample}_semibin_bins"),
        bam_dir= os.path.join(dir_backmapping, "{sample}_cluster_50"),
        sample="{sample}",
    conda:
        os.path.join(dir_env, "semibin_gpu.yaml")
    resources:
        mem_mb =config['resources']['gpujob']['mem_mb'],
        runtime = config['resources']['gpujob']['runtime'],
        partition = "gpu",
        gres = "gpu:1"
    threads:
        config['resources']['gpujob']['threads']
    shell:
        """
        #generating a concatenated fasta file for semibin2 but with only one assembly file 
        #This is being done for renaming the contigs 
        #generating the bins with semibin2
        if [ -f {params.bin_dir}/output_bins/SemiBin_0.fq.gz ]; then
            echo "Already run skipping semibin2 for {params.sample}"
        else
            SemiBin2 single_easy_bin -i {input.assembly} -b {params.bam_dir}/*.bam -o {params.bin_dir} -t {threads} --engine gpu
            touch {output.bins}
        fi
        touch {output.bins}
        """
 
rule semibin_multi_sample_all2all:
    input:
        assembly = os.path.join(dir_assembly,"{sample}.megahit.contigs.fa.gz"),
        bam = os.path.join(dir_backmapping, "{sample}_bam", "done.txt")
    output:
        bins = os.path.join(dir_binning, "{sample}_semibin_bins", "done.txt")
    params:
        bin_dir= os.path.join(dir_binning, "{sample}_semibin_bins"),
        bam_dir= os.path.join(dir_backmapping, "{sample}_bam"),
        sample="{sample}",
    conda:
        os.path.join(dir_env, "semibin_gpu.yaml")
    resources:
        mem_mb =config['resources']['gpujob']['mem_mb'],
        runtime = config['resources']['gpujob']['runtime'],
        partition = "gpu",
        gres = "gpu:1"
    threads:
        config['resources']['long_shortjob']['threads']
    shell:
        """
        #generating a concatenated fasta file for semibin2 but with only one assembly file 
        #This is being done for renaming the contigs 
        #generating the bins with semibin2
        if [ -f {params.bin_dir}/output_bins/SemiBin_0.fq.gz ]; then
            echo "Already run skipping semibin2 for {params.sample}"
        else
            SemiBin2 single_easy_bin -i {input.assembly} -b {params.bam_dir}/*.bam -o {params.bin_dir} -t {threads} --engine gpu
            touch {output.bins}
        fi
        touch {output.bins}
        """

rule merge_semibins:
    input:
        ins = expand(os.path.join(dir_binning, "{sample}_semibin_bins", "done.txt"), sample=sample_names)
    output:
        os.path.join(dir_binning, "all_semibin_bins", "done.txt")
    localrule: True
    params:
        dirs=os.path.join(dir_binning, "all_semibin_bins"),
        bins=os.path.join(dir_binning),
        sample=" ".join(sample_names)
    shell:
        """
        mkdir -p {params.dirs}

        for f in {params.bins}/*_semibin_bins; do
            sample=$(basename "$f" _semibin_bins)

            for file in "$f"/output_bins/*.fa.gz; do
                [ -e "$file" ] || continue

                bn=$(basename "$file")
                newname="${{sample}}_${{bn}}"

                cp "$file" "{params.dirs}/$newname"
            done
        done

        touch {output}
        """