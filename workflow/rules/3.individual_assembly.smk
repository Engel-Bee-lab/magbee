"""
Assembly rules
Illumina paired end reads - individual assembly using megahit
"""
rule megahit_individual_assembly:
    input:
        r1 = os.path.join(dir_hostcleaned,"{sample}_R1.hostcleaned.fastq.gz"),
        r2 = os.path.join(dir_hostcleaned,"{sample}_R2.hostcleaned.fastq.gz")
    output:
        assembly = os.path.join(dir_assembly,"{sample}.megahit.contigs.fa")
    params:
        megahit=os.path.join(dir_assembly,"{sample}.megahit")
    conda:
        os.path.join(dir_env, "megahit.yaml")
    resources:
        mem =config['resources']['bigjob']['mem'],
        time = config['resources']['bigjob']['time']
    threads: 
        config['resources']['bigjob']['cpu']
    shell:
        """
        megahit -1 {input.r1} -2 {input.r2} -o {params.megahit}
        cp {params.megahit}/final.contigs.fa {output.assembly}
        """

rule quast_individual:
    input:
        assembly = os.path.join(dir_assembly,"{sample}.megahit.contigs.fa")
    output:
        report = os.path.join(dir_assembly, "{sample}_quast_output", "report.txt")
    params:
        outdir=os.path.join(dir_assembly, "{sample}_quast_output")
    conda:
        os.path.join(dir_env, "megahit.yaml")
    resources:
        mem =config['resources']['smalljob']['mem'],
        time = config['resources']['smalljob']['time']
    threads: 
        config['resources']['smalljob']['cpu']
    shell:
        """
        quast {input.assembly} -o {params.outdir} --threads {threads} 
        """