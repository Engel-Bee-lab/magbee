"""
Mapping rules
simka rule for picking the 50 assemblies run in the previous workflow
Now backmapping based on the 50 assemblies
"""
import os
import glob 

def read_cluster_map():
    mapping = {}

    for sample in sample_names:
        cluster_file = os.path.join(
            dir_binning,
            f"{sample}_cluster_50/{sample}_backmap_samples.txt"
        )

        with open(cluster_file) as f:
            mapping[sample] = [l.strip() for l in f if l.strip()]

    return mapping

CLUSTER_MAP = read_cluster_map()

rule assembly_index:
    input:
        assembly = os.path.join(dir_assembly,"{sample}.megahit.contigs.fa")
    output:
        index=os.path.join(dir_assembly, "{sample}.mmi"),
    conda:
        os.path.join(dir_env, "minimap2.yaml")
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        minimap2 -d {output.index} {input.assembly}
        """


rule bakckmapping_simka:
    input:
        assembly = os.path.join(dir_assembly,"{sample}.megahit.contigs.fa"),
        index=os.path.join(dir_assembly, "{sample}.mmi"),
        simka_clusters_txt = os.path.join(dir_binning, "{sample}_cluster_50", "{sample}_backmap_samples.txt")
    output:
        bam=os.path.join(dir_binning, "{sample}_cluster_50", "done.txt")
    params:
        csamples=lambda wildcards: CLUSTER_MAP[wildcards.sample],
        reads_path=dir_hostcleaned,
        bam_path=os.path.join(dir_binning, "{sample}_cluster_50"),
        wsample="{sample}"
    conda:
        os.path.join(dir_env, "minimap2.yaml")
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        for sim_sample in {params.csamples}; do
            r1={params.reads_path}/${{sim_sample}}_R1.hostcleaned.fastq.gz
            r2={params.reads_path}/${{sim_sample}}_R2.hostcleaned.fastq.gz

            outbam={params.bam_path}/{params.wsample}_cluster_50_${{sim_sample}}.bam

            minimap2 -ax sr -t {threads} {input.index} $r1 $r2 | samtools view -bS - | samtools sort -o $outbam                
            samtools index $outbam
        done
        touch {output.bam}
        """
