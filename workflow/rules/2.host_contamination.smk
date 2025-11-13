"""
Rules for host contamination removal 
"""
rule host_contamination_removal:
    input:
        r1 = os.path.join(dir_fastp,"{sample}_R1.fastq.gz"),
        r2 = os.path.join(dir_fastp,"{sample}_R2.fastq.gz"),
        host= config['args']['host_seq']
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

        #Mapping to host genome using minimap2 and extracting unmapped reads, versy strict criteria
        # -F 256 to filter secondary alignments, keeps only primary alignments
        # -f 12 to keep reads where both reads are unmapped
    
        minimap2 -ax sr -t {threads} {input.host} {input.r1} {input.r2} \
            | samtools view -b -f 12 -F 256 -@ {threads} -o {params.unmapped_bam} -
        
        samtools fastq -@ {threads} -0 /dev/null -s /dev/null -n \
            -1 >(gzip -c  > {output.r1}) \
            -2 >(gzip -c  > {output.r2}) \
            {params.unmapped_bam}
        
        touch {output.r1}
        touch {output.r2}
        """