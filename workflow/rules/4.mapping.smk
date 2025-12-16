"""
Mapping rules
For now mapping reads to the individual assemblies
all versus all
"""

from glob import glob

rule make_index_of_individual_assembly:
    input:
        assembly = os.path.join(dir_assembly,"{sample}.megahit.contigs.fa")
    output:
        index = os.path.join(dir_binning,"{sample}_index.mmi")
    conda:
        os.path.join(dir_env, "minimap2.yaml")
    resources:
        mem =config['resources']['smalljob']['mem'],
        time = config['resources']['smalljob']['time']
    threads: 
        config['resources']['smalljob']['cpu']
    shell:
        """
        if [ -d {output.index} ]; then
            echo "Index already exists. Skipping..."
        else 
            minimap2 {input.assembly} -d {output.index}
        fi
        """

rule map_reads_to_individual_assembly:
    input:
        index = os.path.join(dir_binning, "{sample}_index.mmi"),
        assemble = os.path.join(dir_assembly,"{sample}.megahit.contigs.fa")
    params:
        bins= os.path.join(dir_binning),
        samples= os.path.join(dir_binning, "{sample}_bam"),
        s=expand("{sample}", sample=sample_names),
        ids="{sample}",
        reads=os.path.join(dir_hostcleaned)
    output:
        txt = os.path.join(dir_binning, "{sample}_bam", "done.txt")
    conda:
        os.path.join(dir_env, "minimap2.yaml")
    resources:
        mem =config['resources']['smalljob']['mem'],
        time = config['resources']['smalljob']['time']
    threads: 
        config['resources']['smalljob']['cpu']
    shell:
        """
        if [ -d {output.txt} ]; then
            echo "Mapping already done. Skipping..."
        else
            mkdir -p {params.samples}
            for f in {params.s}; do
                minimap2 -ax sr -t {threads} {input.index} {params.reads}/"$f"_R1.hostcleaned.fastq.gz {params.reads}/"$f"_R2.hostcleaned.fastq.gz | samtools view -bS - | samtools sort -o {params.samples}/{params.ids}Ref_"$f"Reads.bam
                samtools index {params.samples}/{params.ids}Ref_"$f"Reads.bam
            done
            touch {output.txt}
        fi
        """