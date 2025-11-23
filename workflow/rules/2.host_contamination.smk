"""
Rules for host contamination removal 
"""
rule host_mapping:
    input:
        r1 = os.path.join(dir_fastp,"{sample}_R1.fastq.gz"),
        r2 = os.path.join(dir_fastp,"{sample}_R2.fastq.gz"),
        host= config['args']['host_seq']
    output:
        r1 = os.path.join(dir_hostcleaned,"{sample}_R1.hostcleaned.fastq.gz"),
        r2 = os.path.join(dir_hostcleaned,"{sample}_R2.hostcleaned.fastq.gz"),
        mr1= os.path.join(dir_hostcleaned,"{sample}_R1.mapped.fastq.gz"),
        mr2= os.path.join(dir_hostcleaned,"{sample}_R2.mapped.fastq.gz")
    params:
        all_bam=temporary(os.path.join(dir_hostcleaned,"{sample}_temp.bam")),
        unmapped_bam = os.path.join(dir_hostcleaned,"{sample}_unmapped.bam"),
        mapped_bam = os.path.join(dir_hostcleaned,"{sample}_mapped.bam")
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
        minimap2 -ax sr -t {threads} {input.host} {input.r1} {input.r2} \
            | samtools view -b -@ {threads} -o {params.all_bam} -

        samtools view -b -F 4 -f 12 -@ {threads} {params.all_bam} -o {params.unmapped_bam}
        samtools view -b -F 4 -F 256 -@ {threads} {params.all_bam} -o {params.mapped_bam} 

        #unmapped reads to fastq
        samtools fastq -@ {threads} -0 /dev/null -s /dev/null -n \
            -1 >(gzip -c  > {output.r1}) \
            -2 >(gzip -c  > {output.r2}) \
            {params.unmapped_bam}
        
        #mapped reads to fastq (if needed)
        samtools fastq -@ {threads} -0 /dev/null -s /dev/null -n \
            -1 >(gzip -c  > {output.mr1}) \
            -2 >(gzip -c  > {output.mr2}) \
            {params.mapped_bam}

        touch {output.r1}
        touch {output.r2}
        touch {output.mr1}
        touch {output.mr2}
        """
