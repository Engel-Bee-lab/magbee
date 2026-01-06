"""
Extra steps with the host reads 
    - map to the the reference mitogenome
"""

rule host_mapping:
    input:
        mr1 = os.path.join(dir_hostcleaned,"{sample}_R1.mapped.fastq.gz"),
        mr2 = os.path.join(dir_hostcleaned,"{sample}_R2.mapped.fastq.gz"),
        host= config['extra_db']['mitogenome']
    output:
        mr1_mt = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mt_R1.mapped.fastq.gz"),
        mr2_mt = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mt_R2.mapped.fastq.gz"),
    params:
        mapped_bam=os.path.join(dir_hostcleaned,"{sample}_mapped.bam"),
        stats=os.path.join(dir_hostcleaned,"{sample}_bamstats.txt")
    conda:
        os.path.join(dir_env, "minimap2.yaml")
    resources:
        mem =config['resources']['smalljob']['mem'],
        time = config['resources']['smalljob']['time']
    threads: 
        config['resources']['smalljob']['cpu']
    shell:
        """
        set -euo pipefail
        if [ -f {output.mr1_mt} ] && [ -f {output.mr2_mt} ]; then
            echo "Host mitogenome mapping already done. Skipping..."
            exit 0
        else
            minimap2 -ax sr -t {threads} {input.host} {input.mr1} {input.mr2} \
                | samtools view -b -F 4 -@ {threads} -o {params.mapped_bam} -
            samtools flagstat {params.mapped_bam} > {params.stats}

            samtools fastq -@ {threads} -0 /dev/null -s /dev/null -n \
                -1 >(gzip -c  > {output.mr1_mt}) \
                -2 >(gzip -c  > {output.mr2_mt}) \
                {params.mapped_bam}
            
            touch {output.mr1_mt}
            touch {output.mr2_mt}
        
        fi
        """