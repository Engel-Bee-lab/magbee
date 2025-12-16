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
        mem =config['resources']['assemblyjob']['mem'],
        time = config['resources']['assemblyjob']['time']
    threads: 
        config['resources']['assemblyjob']['cpu']
    shell:
        """
        if [ -d {params.megahit} ]; then
            echo "Megahit already run."
            if [ ! -f {params.megahit}/final.contigs.fa ]; then
                echo "But final contigs file not found, rerunning megahit."
                megahit -1 {input.r1} -2 {input.r2} -o {params.megahit} -t {threads} -m {resources.mem} -f
            else
                echo "Final contigs file found."
                cp {params.megahit}/final.contigs.fa {output.assembly}
            fi
        else
            megahit -1 {input.r1} -2 {input.r2} -o {params.megahit} -t {threads} -m {resources.mem}
            cp {params.megahit}/final.contigs.fa {output.assembly}
        fi
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