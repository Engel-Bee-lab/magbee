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
        mem_mb =config['resources']['assemblyjob']['mem_mb'],
        runtime = config['resources']['assemblyjob']['runtime']
    threads: 
        config['resources']['assemblyjob']['threads']
    shell:
        """
        if [ -d {params.megahit} ]; then
            echo "Megahit already run."
            if [ ! -f {params.megahit}/final.contigs.fa ]; then
                echo "But final contigs file not found, rerunning megahit."
                megahit -1 {input.r1} -2 {input.r2} -o {params.megahit} -t {threads} -m {resources.mem_mb} -f
            else
                echo "Final contigs file found."
                cp {params.megahit}/final.contigs.fa {output.assembly}
            fi
        else
            megahit -1 {input.r1} -2 {input.r2} -o {params.megahit} -t {threads} -m {resources.mem_mb}
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
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        quast {input.assembly} -o {params.outdir} --threads {threads}
        """