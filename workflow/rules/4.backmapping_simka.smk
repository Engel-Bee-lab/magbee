"""
Mapping rules
simka rule for picking the 50 assemblies run in the previous workflow
Now backmapping based on the 50 assemblies
"""
import os
import glob 

############################################
# Helper: read cluster file SAFELY (lazy)
############################################
def get_cluster_samples(wildcards, input):
    with open(input.simka_clusters_txt) as f:
        return [l.strip() for l in f if l.strip()]

  
#CLUSTER_MAP = read_cluster_map()

############################################
# Rule: index assembly
############################################
rule assembly_index:
    input:
        assembly = os.path.join(dir_assembly, "{sample}.megahit.contigs.fa.gz")
    output:
        index = os.path.join(dir_assembly, "{sample}.mmi")
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

############################################
# Rule: backmapping reads to assembly
############################################
rule backmapping_simka:
    input:
        assembly = os.path.join(dir_assembly, "{sample}.megahit.contigs.fa.gz"),
        index = os.path.join(dir_assembly, "{sample}.mmi"),
        simka_clusters_txt = os.path.join(
            dir_binning,
            "{sample}_cluster_50",
            "{sample}_backmap_samples.txt"
        )
    output:
        done = os.path.join(dir_backmapping, "{sample}_cluster_50", "done.txt")
    params:
        csamples = get_cluster_samples,
        reads_path = dir_hostcleaned,
        bam_path = os.path.join(dir_backmapping, "{sample}_cluster_50"),
        wsample = "{sample}"
    conda:
        os.path.join(dir_env, "minimap2.yaml")
    threads:
        config['resources']['long_shortjob']['threads']
    resources:
        mem_mb = config['resources']['long_shortjob']['mem_mb'],
        runtime = config['resources']['long_shortjob']['runtime']
    shell:
        """
        mkdir -p {params.bam_path}

        for sim_sample in {params.csamples}; do
            r1={params.reads_path}/${{sim_sample}}_R1.hostcleaned.fastq.gz
            r2={params.reads_path}/${{sim_sample}}_R2.hostcleaned.fastq.gz

            outbam={params.bam_path}/{params.wsample}_cluster_50_${{sim_sample}}.bam

            minimap2 -ax sr -t {threads} {input.index} $r1 $r2 \
                | samtools view -bS - \
                | samtools sort -o $outbam

            samtools index $outbam
        done

        touch {output.done}
        """