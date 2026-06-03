"""
Rules to concatenate the assemblies and then map all the reads to this assembly
"""
import os
import glob 

rule rename_contigs:
    input:
        os.path.join(dir_assembly, "{sample}.megahit.contigs.fa.gz"), sample=sample_names
    output:
        os.path.join(dir_temp, "{sample}.megahit.contigs.fa.gz"), sample=sample_names
    params:
        sample=lambda wc: wc.sample
    localrule: True
    shell:
        """
        zcat {input} | \
        awk -v prefix="{wildcards.sample}" '
            /^>/ {{print ">"prefix"_"substr($0,2); next}}
            {{print}}
        ' | gzip > {output}
        """

rule concatenate_assemblies:
    input:
        expand(os.path.join(dir_temp, "{sample}.megahit.contigs.fa.gz"), sample=sample_names)
    output:
       concatenated_assembly = os.path.join(dir_temp, "concatenated_assemblies.fa.gz")
    localrule: True
    shell:
        """
        zcat {input} | gzip > {output}
        """

rule assembly_index:
    input:
        assembly = os.path.join(dir_temp, "concatenated_assemblies.fa.gz")
    output:
        index = os.path.join(dir_temp, "concatenated_assemblies.mmi")
    conda:
        os.path.join(dir_env, "minimap2.yaml")
    threads:
        config['resources']['bigjob']['threads']
    resources:
        mem_mb = config['resources']['bigjob']['mem_mb'],
        runtime = config['resources']['bigjob']['runtime']
    shell:
        """
        minimap2 -d {output.index} {input.assembly}
        """

rule map_reads_to_assembly:
    input:
        assembly = os.path.join(dir_temp, "concatenated_assemblies.fa.gz"),
        index = os.path.join(dir_temp, "concatenated_assemblies.mmi"),
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
        