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
        config['resources']['smalljob']['threads']
    resources:
        mem_mb = config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    shell:
        """
        minimap2 -d {output.index} {input.assembly}
        """

rule map_reads_to_assembly:
    input:
        assembly = os.path.join(dir_assembly, "concatenated_assemblies.fa.gz"),
        index = os.path.join(dir_assembly, "concatenated_assemblies.mmi"),
    output:
        txt = os.path.join(dir_backmapping, "{sample}_bam", "done.txt")
    params:
        samples="{sample}",
        reads=os.path.join(dir_hostcleaned)
    conda:
        os.path.join(dir_env, "minimap2.yaml")
    threads:
        config['resources']['smalljob']['threads']
    resources:
        mem_mb = config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    shell:
        """
        if [ -d {output.txt} ]; then
            echo "Mapping already done. Skipping..."
        else
            mkdir -p {params.samples}
            minimap2 -ax sr -t {threads} {input.index} {params.reads}/{params.samples}_R1.hostcleaned.fastq.gz {params.reads}/{params.samples}_R2.hostcleaned.fastq.gz \
                | samtools view -bS - | samtools sort -o {params.samples}/{params.samples}.bam
                samtools index {params.samples}/{params.samples}.bam
            done
            touch {output.txt}
        fi
        """
        