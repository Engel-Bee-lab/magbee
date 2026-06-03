"""
Rules to concatenate the assemblies and then map all the reads to this assembly
"""
import os
import glob 

rule concatenate_assemblies:
    input:
        expand(os.path.join(dir_assembly, "{sample}.megahit.contigs.fa.gz"), sample=sample_names)
    output:
       concatenated_assembly = os.path.join(dir_assembly, "concatenated_assemblies.fa.gz")
    localrule: True
    shell:
        """
        cat {input} > {output}
        """

rule assembly_index:
    input:
        assembly = os.path.join(dir_assembly, "concatenated_assemblies.fa.gz")
    output:
        index = os.path.join(dir_assembly, "concatenated_assemblies.mmi")
    conda:
        os.path.join(dir_env, "minimap2.yaml")
    threads:
        config['resources']['bigjob']['threads']
    resources:
        mem_mb = config['resources']['bigjob']['mem_mb'],
        runtime = config['resources']['bigjob']['runtime']
    shell:
        """
        minimap2 {input.assembly} -d {output.index}
        """

rule map_reads_to_assembly:
    input:
        assembly = os.path.join(dir_assembly, "concatenated_assemblies.fa.gz"),
        index = os.path.join(dir_assembly, "concatenated_assemblies.mmi"),
    output:
        txt = os.path.join(dir_backmapping, "{sample}_bam", "done.txt")
    params:
        r1=lambda wc: os.path.join(dir_hostcleaned, f"{wc.sample}_R1.hostcleaned.fastq.gz"),
        r2=lambda wc: os.path.join(dir_hostcleaned, f"{wc.sample}_R2.hostcleaned.fastq.gz"),
        outdir=lambda wc: os.path.join(dir_backmapping, f"{wc.sample}_bam")
    conda:
        os.path.join(dir_env, "minimap2.yaml")
    threads:
        config['resources']['bigjob']['threads']
    resources:
        mem_mb = config['resources']['bigjob']['mem_mb'],
        runtime = config['resources']['bigjob']['runtime']
    shell:
        """
        mkdir -p {params.outdir}

        minimap2 -ax sr -t {threads} {input.index} \
            {params.r1} {params.r2} \
        | samtools view -bS - \
        | samtools sort -o {params.outdir}/{wildcards.sample}.bam

        samtools index {params.outdir}/{wildcards.sample}.bam

        touch {output.txt}
        """
        