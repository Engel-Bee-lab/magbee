"""
Assembly rules
Illumina paired end reads - co-assembly using megahit
"""
rule megahit_coassembly:
    input:
        r1 = expand(os.path.join(dir_hostcleaned,"{sample}_R1.hostcleaned.fastq.gz"), sample=samples_names),
        r2 = expand(os.path.join(dir_hostcleaned,"{sample}_R2.hostcleaned.fastq.gz"), sample=samples_names)
    output:
        assembly = os.path.join(dir_assembly,"coassembly.megahit.contigs.fa")
    conda:
        os.path.join(dir_env, "megahit.yaml")
    resources:
        mem_mb =config['resources']['bigjob']['mem_mb'],
        runtime = config['resources']['bigjob']['runtime']
    threads: 
        config['resources']['bigjob']['threads']
    shell:
        """
        cat {input.r1} > {dir_assembly}/coassembly_R1.fastq.gz
        cat {input.r2} > {dir_assembly}/coassembly_R2.fastq.gz

        megahit -1 {dir_assembly}/coassembly_R1.fastq.gz -2 {dir_assembly}/coassembly_R2.fastq.gz -o {dir_assembly}/coassembly.megahit --min-contig-len 1000
        cp {dir_assembly}/coassembly.megahit/final.contigs.fa {output.assembly}
        """

rule quast:
    input:
        assembly = os.path.join(dir_assembly,"coassembly.megahit.contigs.fa")
    output:
        report = os.path.join(dir_assembly,"quast_report.txt")
    conda:
        os.path.join(dir_env, "megahit.yaml")
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        quast {input.assembly} -o {dir_assembly}/quast_output --threads {threads} > {output.report}
        """
