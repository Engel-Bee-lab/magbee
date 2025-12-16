"""
Rules for host contamination removal 
"""
rule host_mapping:
    input:
        r1 = os.path.join(dir_fastp,"{sample}_R1.fastq.gz"),
        r2 = os.path.join(dir_fastp,"{sample}_R2.fastq.gz"),
        host= config['args']['host_seq']
    output:
        all_bam=os.path.join(dir_hostcleaned,"{sample}_temp.bam"),
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
        if [ -f {output.all_bam} ] && [ -f {output.stats} ]; then
            echo "Host mapping already done. Skipping..."
            exit 0
        else
            minimap2 -ax sr -t {threads} {input.host} {input.r1} {input.r2} \
                | samtools view -b -@ {threads} -o {output.all_bam} -
            samtools flagstat {output.all_bam} > {output.stats}
        fi
        """

"""
Rules for host contamination removal 
"""
rule unmapped_reads:
    input:
        r1 = os.path.join(dir_fastp,"{sample}_R1.fastq.gz"),
        r2 = os.path.join(dir_fastp,"{sample}_R2.fastq.gz"),
        all_bam=os.path.join(dir_hostcleaned,"{sample}_temp.bam")
    output:
        r1 = os.path.join(dir_hostcleaned,"{sample}_R1.hostcleaned.fastq.gz"),
        r2 = os.path.join(dir_hostcleaned,"{sample}_R2.hostcleaned.fastq.gz"),
    params:
        unmapped_bam = os.path.join(dir_hostcleaned,"{sample}_unmapped.bam")
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
        if [ -f {output.r1} ] && [ -f {output.r2} ]; then
            echo "Unmapped reads already exist. Skipping..."
            exit 0
        else
            # get unmapped reads
            samtools view -b -f 4 -@ {threads} -o {params.unmapped_bam} {input.all_bam}

            #unmapped reads to fastq
            samtools fastq -@ {threads} -0 /dev/null -s /dev/null -n \
                -1 >(gzip -c  > {output.r1}) \
                -2 >(gzip -c  > {output.r2}) \
                {params.unmapped_bam}

            touch {output.r1}
            touch {output.r2}
        fi
        """

"""
Rules for host reads extraction 
"""
rule host_mapped_reads:
    input:
        r1 = os.path.join(dir_fastp,"{sample}_R1.fastq.gz"),
        r2 = os.path.join(dir_fastp,"{sample}_R2.fastq.gz"),
        all_bam=os.path.join(dir_hostcleaned,"{sample}_temp.bam")
    output:
        mr1 = os.path.join(dir_hostcleaned,"{sample}_R1.mapped.fastq.gz"),
        mr2 = os.path.join(dir_hostcleaned,"{sample}_R2.mapped.fastq.gz"),
    params:
        mapped_bam = os.path.join(dir_hostcleaned,"{sample}_mapped.bam"),
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
        if [ -f {output.mr1} ] && [ -f {output.mr2} ]; then
            echo "Mapped reads already exist. Skipping..."
            exit 0
        else
            samtools view -b -F 4 -@ {threads} -o {params.mapped_bam} {input.all_bam}

            #unmapped reads to fastq
            samtools fastq -@ {threads} -0 /dev/null -s /dev/null -n \
                -1 >(gzip -c  > {output.mr1}) \
                -2 >(gzip -c  > {output.mr2}) \
                {params.mapped_bam}
            
            touch {output.mr1}
            touch {output.mr2}
        fi
        """