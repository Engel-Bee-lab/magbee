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
        os.path.join(dir_env, "metaspades.yaml")
    resources:
        mem =config['resources']['largejob']['mem'],
        time = config['resources']['largejob']['time']
    threads: 
        config['resources']['largejob']['cpu']
    shell:
        """
        spades.py --meta -1 {input.r1} -2 {input.r2} -o {params.megahit} -t {threads} -m {resources.mem}
        cp {params.megahit}/contigs.fasta {output.assembly}
        """